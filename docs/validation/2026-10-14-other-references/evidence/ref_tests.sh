#!/usr/bin/env bash
# bonefish vs minimap2 on other references, genome-wide HG001 sample (161,762 reads):
#   T2T-CHM13v2.0 (UCSC hs1), and GRCh38 with ALT contigs and decoys (NCBI full_plus_hs38d1 analysis set)
#   with the ALT list given to both aligners (--alt). Index builds and mapping timed; mapping statistics,
#   primaries on ALT contigs, and MAPQ of primaries in the MHC (chr6:28.5-33.4 Mb).
set -uo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
OUT="$ROOT/results/reftests"; mkdir -p "$OUT/runs" "$OUT/stats"
X="$ROOT/data/ref/extra"; IDX="$ROOT/data/index/extra"; mkdir -p "$IDX"
READS="$ROOT/data/reads/hg001_wg_1x.fq.gz"; T=16
BF="$ROOT/bonefish/target/release/bonefish"; MM2=$HOME/micromamba/envs/ontaln/bin/minimap2
TV="/usr/bin/time -v -o"
log() { echo "[$(date +%T)] $*" | tee -a "$OUT/run.log"; }
# wait for the wf-human-variation runs (CPU) to finish
until grep -q ALLDONE "$ROOT/results/wfhv/run.log" 2>/dev/null && ! pgrep -f "nextflow.*wf-human-variation" > /dev/null; do sleep 60; done
log "start"
[[ -s $X/hs1.fa ]] || zcat $X/hs1.fa.gz > $X/hs1.fa
[[ -s $X/GRCh38_alt.fa ]] || zcat $X/GRCh38_full_plus_hs38d1.fna.gz | awk '/^>/{print $1; next} {print}' > $X/GRCh38_alt.fa
grep '^>' $X/GRCh38_alt.fa | sed 's/>//' | grep '_alt$' > $X/GRCh38_alt.list
log "ALT contigs: $(wc -l < $X/GRCh38_alt.list)"
for r in hs1 GRCh38_alt; do
  [[ -s $IDX/$r.bonefish.raw ]] || { log "index bonefish $r"; $TV $OUT/runs/idx_bonefish_$r.time $BF -t $T -d $IDX/$r.bonefish.raw --idx-raw $X/$r.fa 2>/dev/null; }
  [[ -s $IDX/$r.mm2.mmi ]] || { log "index minimap2 $r"; $TV $OUT/runs/idx_minimap2_$r.time $MM2 -t $T -x lr:hq -d $IDX/$r.mm2.mmi $X/$r.fa 2>/dev/null; }
done
run() { # name, command
  local n=$1; shift
  [[ -s $OUT/stats/$n/bam_summary.tsv ]] && return
  log "map $n"
  bash -c "$TV $OUT/runs/$n.time $* 2>$OUT/runs/$n.err | cut -f1-6,12- | grep -v '^@' > $OUT/runs/$n.cols"
  mkdir -p $OUT/stats/$n; python3 $ROOT/bench/bam_summary.py $OUT/stats/$n < $OUT/runs/$n.cols
  log "done $n $(grep -E 'Elapsed' $OUT/runs/$n.time | sed 's/.*: //')"
}
cat $IDX/hs1.bonefish.raw > /dev/null; run chm13_bonefish "$BF -t $T -a $IDX/hs1.bonefish.raw $READS"
cat $IDX/hs1.mm2.mmi > /dev/null; run chm13_minimap2 "$MM2 -t $T -ax lr:hq $IDX/hs1.mm2.mmi $READS"
cat $IDX/GRCh38_alt.bonefish.raw > /dev/null; run alt_bonefish "$BF -t $T -a --alt $X/GRCh38_alt.list $IDX/GRCh38_alt.bonefish.raw $READS"
cat $IDX/GRCh38_alt.mm2.mmi > /dev/null; run alt_minimap2 "$MM2 -t $T -ax lr:hq --alt $X/GRCh38_alt.list $IDX/GRCh38_alt.mm2.mmi $READS"
# ALT-specific statistics from the kept columns
for n in alt_bonefish alt_minimap2; do
  awk -F'\t' -v n=$n 'int($2/256)%2==0 && int($2/4)%2==0 && int($2/2048)%2==0 { p++; if ($3 ~ /_alt$/) a++; if ($3 ~ /decoy/) d++; if ($3=="chr6" && $4>=28500000 && $4<=33400000) { m++; if ($5>=20) m20++ } }
    END { printf "%s\tprimaries=%d\ton_alt=%d\ton_decoy=%d\tMHC_primaries=%d\tMHC_MAPQ>=20=%d (%.2f%%)\n", n, p, a, d, m, m20, 100*m20/(m?m:1) }' $OUT/runs/$n.cols
done | tee $OUT/alt_stats.tsv
log ALLDONE
