#!/usr/bin/env bash
# CpG methylation from bonefish vs minimap2 BAMs: HG002 chr20 reads (PAW70337, dorado sup v5.0.0 with
# 5mC/5hmC and 6mA calls) aligned genome-wide by both, then modkit pileup on chr20 (--cpg
# --combine-strands) and a per-site comparison (bench/modkit_compare.py).
# Usage: bench/modkit_compare.sh [OUT_DIR]
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
OUT="${1:-$ROOT/results/modkit}"; mkdir -p "$OUT"
T=16
REF="$ROOT/data/ref/GRCh38_no_alt.fa"
SRC="$ROOT/data/reads/PAW70337.sparse.bam"
SAMTOOLS=$HOME/micromamba/envs/vareval/bin/samtools
MM2=$HOME/micromamba/envs/ontaln/bin/minimap2
MODKIT=$HOME/micromamba/envs/alnbench/bin/modkit
BF="$ROOT/bonefish/target/release/bonefish"
IDX="$ROOT/data/index/bench"
log() { echo "[$(date +%T)] $*" | tee -a "$OUT/run.log"; }
TV="/usr/bin/time -v -o"

if [[ ! -s $OUT/reads.bam ]]; then
  log "extract chr20 primaries"
  $SAMTOOLS view -b -F 0x900 -o $OUT/reads.bam $SRC chr20
fi
if [[ ! -s $OUT/reads.fq ]]; then
  $SAMTOOLS fastq -@4 -T MM,ML,MN $OUT/reads.bam > $OUT/reads.fq 2>/dev/null
fi
log "reads: $(awk 'NR%4==2{n++; b+=length($0)} END{print n" reads, "b" bases"}' $OUT/reads.fq)"

if [[ ! -s $OUT/minimap2.bam.csi && ! -s $OUT/minimap2.bam.bai ]]; then
  log "minimap2"
  $TV $OUT/minimap2.time bash -c "$MM2 -t $T -ax lr:hq -y $IDX/mm2.lrhq.mmi $OUT/reads.fq 2>$OUT/minimap2.err | $SAMTOOLS sort -@8 -m 1G --write-index -o $OUT/minimap2.bam -"
fi
if [[ ! -s $OUT/bonefish.bam.bai ]]; then
  log "bonefish (FASTQ with tags)"
  $TV $OUT/bonefish.time $BF -t $T -y -o $OUT/bonefish.bam $IDX/bonefish.raw $OUT/reads.fq 2>$OUT/bonefish.err
fi
if [[ ! -s $OUT/bonefish_bamin.bam.bai ]]; then
  log "bonefish (BAM input)"
  $TV $OUT/bonefish_bamin.time $BF -t $T -y -o $OUT/bonefish_bamin.bam $IDX/bonefish.raw $OUT/reads.bam 2>$OUT/bonefish_bamin.err
fi
for a in minimap2 bonefish bonefish_bamin; do
  if [[ ! -s $OUT/$a.cpg.bed ]]; then
    log "modkit pileup $a"
    $TV $OUT/$a.modkit.time $MODKIT pileup $OUT/$a.bam $OUT/$a.cpg.bed --modified-bases 5mC 5hmC --cpg --ref $REF --combine-strands --region chr20 \
      -t $T --log-filepath $OUT/$a.modkit.log 2>$OUT/$a.modkit.err
  fi
done
log ALLDONE
