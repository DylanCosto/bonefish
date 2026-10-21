#!/usr/bin/env bash
# RNA benchmark for one bonefish build (or minimap2 with BIN=minimap2): speed on 200k real reads of HG002
# direct RNA (RNA004 sup, -x rna) and cDNA (PCS114 hac, -x cdna), stage timers, agreement with
# minimap2 2.31 (records, intron chains), GENCODE junction share, and simulated-read truth.
# Usage: rna_bench.sh LABEL BIN [extra bonefish args...]; results in results/rna/runs/LABEL, one line
# appended to results/rna/ledger.tsv. Env: T (threads, 16), SETS (default "drna200k cdna200k"), NOSIM=1.
set -uo pipefail
L=$1; BIN=$2; shift 2; EXTRA="$*"
R=${BENCH}; D=$R/data/rna; O=$R/results/rna/runs/$L; M=$R/results/rna/mm2
export PATH=$HOME/micromamba/envs/vareval/bin:$HOME/micromamba/envs/ontaln/bin:$PATH
T=${T:-16}; SETS=${SETS-"drna200k cdna200k"}
mkdir -p $O; cd $O
declare -A READS=([drna200k]=$D/drna_hg002_200k.fq.gz [cdna200k]=$D/cdna_hg002_200k.fq.gz
                  [sim_drna]=$D/sim/sim_drna.fq.gz [sim_cdna]=$D/sim/sim_cdna.fq.gz
                  [drna]=$D/drna_hg002.fq.gz [cdna]=$D/cdna_hg002.fq.gz)
# The indexes do not fit in the page cache together: drop all but the one about to be used.
evict_others() {
  python3 -c "
import os, sys
for f in sys.argv[2:]:
    if f != sys.argv[1] and os.path.exists(f):
        fd = os.open(f, os.O_RDONLY); os.posix_fadvise(fd, 0, 0, os.POSIX_FADV_DONTNEED); os.close(fd)
" "$1" $D/bf.rna.raw $D/bf.cdna.raw $D/mm2.splice.mmi $D/mm2.splice_k14.mmi
}
preset() { [[ $1 == *drna* ]] && echo rna || echo cdna; }
map() { # set -> $set.sam, $set.time, $set.err
  local s=$1 p; p=$(preset $s)
  local y=""; [[ $s != sim_* ]] && y="-y"
  if [[ $BIN == minimap2 ]]; then
    local idx=$D/mm2.splice.mmi k=""; [[ $p == rna ]] && { idx=$D/mm2.splice_k14.mmi; k="-uf -k14"; }
    evict_others $idx; cat $idx > /dev/null
    /usr/bin/time -f "%e %U %S %M" -o $s.time minimap2 -ax splice $k $y -t $T $idx ${READS[$s]} > $s.sam 2> $s.err
  else
    evict_others $D/bf.$p.raw; cat $D/bf.$p.raw > /dev/null   # page this index in, the others out
    /usr/bin/time -f "%e %U %S %M" -o $s.time $BIN -x $p -a $y -v -t $T $EXTRA $D/bf.$p.raw ${READS[$s]} > $s.sam 2> $s.err
  fi
}
line="$L"
for s in $SETS; do
  map $s
  read -r w u sy m < $s.time
  cpu=$(echo "$u + $sy" | bc)
  line+="\t$s wall=$w cpu=$cpu rss=$((m / 1024))MB"
  if [[ -s $M/${s}.bam ]]; then
    python3 $R/bench/compare_primary.py $M/$s.bam $s.sam > $s.compare.txt
    python3 $R/bench/intron_chains.py $D/gencode.v49.bed12 $M/$s.bam $s.sam > $s.chains.txt
    ident=$(awk '/identical:/{print $3}' $s.compare.txt | tr -d '()')
    chain=$(grep -o 'identical intron chain [0-9,]* ([0-9.]*%)' $s.chains.txt | grep -o '[0-9.]*%')
    ann=$(sed -n 2p $s.chains.txt | grep -o 'annotated [0-9.]*%' | grep -o '[0-9.]*%')
    line+=" ident=$ident chain=$chain annot=$ann"
  fi
done
if [[ -z ${NOSIM:-} ]]; then
  for s in sim_cdna sim_drna; do
    map $s
    python3 $R/bench/rna_sim_truth.py $D/sim/tx.truth.tsv ${READS[$s]} $s.sam > $s.truth.txt
    line+="\t$s $(tail -1 $s.truth.txt | cut -d: -f2- | sed 's/  */ /g')"
  done
fi
rm -f $O/*.sam.tmp
echo -e "$line" | tee -a $R/results/rna/ledger.tsv
