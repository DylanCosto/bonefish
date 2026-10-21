#!/usr/bin/env bash
# Metagenomics (Zymo fecal reference vs its 8-flowcell metaMDBG assembly) and duplex (HG002 vs GRCh38)
# checks: bonefish against minimap2 2.31 lr:hq on the same reads. Run after meta_duplex_fetch.sh.
set -uo pipefail
R=${BENCH}; B=$R/bench; O=$R/results/meta_duplex; mkdir -p $O; cd $O
export PATH=$HOME/micromamba/envs/vareval/bin:$HOME/micromamba/envs/ontaln/bin:$PATH
BF=${BF:-$R/bin/bonefish-v017}; T=16
echo "bonefish: $($BF --version); minimap2: $(minimap2 --version)" > versions.txt
run() { # name idx_bf idx_mm reads
  local n=$1
  [[ -s $n.bf.bam ]] || /usr/bin/time -f "%e s wall, %U+%S CPU-s, %M KB" -o $n.bf.time \
      $BF -a -y -t $T $2 $4 2> $n.bf.err | samtools view -b -o $n.bf.bam -
  [[ -s $n.mm.bam ]] || /usr/bin/time -f "%e s wall, %U+%S CPU-s, %M KB" -o $n.mm.time \
      minimap2 -ax lr:hq -y -t $T $3 $4 2> $n.mm.err | samtools view -b -o $n.mm.bam -
}
M=$R/data/meta
run meta $M/meta.bonefish.raw $M/meta.mm2.mmi $M/fecal_reads.fq.gz
D=$R/data/duplex/11_15_22_R1041_Duplex_HG002_3_Dorado_v0.1.1_400bps_sup_stereo_duplex_pass.fastq.gz
run duplex $R/data/index/bonefish.lrhq.raw $R/data/index/mm2.lrhq.mmi $D
for n in meta duplex; do
  circ=""; [[ $n == meta ]] && circ=$M/bgzip_contigs.fa.gz.circular.list
  for a in bf mm; do
    echo "== $n $a ($(cat $n.$a.time))"
    python3 $B/aln_stats.py $n.$a.bam $n.$a.counts.tsv $circ
  done > $n.stats.txt
  python3 $B/compare_primary.py $n.mm.bam $n.bf.bam > $n.compare.txt
done
python3 $B/abundance_compare.py meta.mm.counts.tsv meta.bf.counts.tsv $M/das_DASTool_contig2bin.tsv > meta.abundance.txt
echo DONE > done.txt
