#!/usr/bin/env bash
# Transcriptome-mapping benchmark: one build (BIN=minimap2 for the -x map-ont -N 100 --eqx baseline)
# on HG002 cDNA / direct RNA 200k (speed, oarfish counts vs minimap2's) and simulated reads (truth).
# Usage: tx_bench.sh LABEL BIN [extra args]; results in results/tx_quant/runs/LABEL, line in ledger.tsv.
set -uo pipefail
L=$1; BIN=$2; shift 2; EXTRA="$*"
R=${BENCH}; D=$R/data/rna; X=$D/tx; O=$R/results/tx_quant/runs/$L; B=$R/results/tx_quant/runs/mm2
export PATH=$HOME/micromamba/envs/vareval/bin:$HOME/micromamba/envs/ontaln/bin:$HOME/micromamba/envs/quant/bin:$PATH
T=16; mkdir -p $O; cd $O
declare -A READS=([cdna]=$D/cdna_hg002_200k.fq.gz [drna]=$D/drna_hg002_200k.fq.gz [sim_cdna]=$D/sim/sim_cdna.fq.gz [sim_drna]=$D/sim/sim_drna.fq.gz)
if [[ $BIN == minimap2 ]]; then IDX=$X/tx.mm2.mmi; else IDX=${IDX:-$X/tx.tx.raw}; fi
cat $IDX > /dev/null
line="$L"
for s in cdna drna sim_cdna sim_drna; do
  if [[ $BIN == minimap2 ]]; then
    /usr/bin/time -f "%e %U %S %M" -o $s.time minimap2 -ax map-ont -N 100 --eqx -t $T $IDX ${READS[$s]} 2> $s.err | samtools view -b -o $s.bam -
  else
    /usr/bin/time -f "%e %U %S %M" -o $s.time $BIN -a --eqx -t $T $EXTRA $IDX ${READS[$s]} 2> $s.err | samtools view -b -o $s.bam -
  fi
  read -r w u sy m < $s.time
  line+="\t$s wall=$w cpu=$(echo "$u + $sy" | bc)"
  if [[ $s == sim_* ]]; then
    line+=" $(python3 $R/bench/tx_sim_truth.py ${READS[$s]} $s.bam | cut -d: -f2- | sed 's/  */ /g')"
  else
    # oarfish accepts only BAMs whose @PG names a validated aligner; relabel bonefish as rammap (its upstream).
    samtools view -H $s.bam | sed 's/^@PG\tID:bonefish\tPN:bonefish/@PG\tID:rammap\tPN:rammap/' > hdr.sam
    samtools reheader hdr.sam $s.bam > $s.q.bam && rm hdr.sam
    oarfish -j $T -a $s.q.bam -o $s.oarfish --filter-group nanocount-filters > $s.oarfish.log 2>&1; rm -f $s.q.bam
    line+=" unmapped=$(samtools view -c -f 4 $s.bam)"
    [[ -s $B/$s.oarfish.quant && $L != mm2 ]] && line+=" $(python3 $R/bench/tx_quant_compare.py $B/$s.oarfish.quant $s.oarfish.quant | tr '\n' ' ' | sed 's/  */ /g')"
  fi
done
echo -e "$line" | tee -a $R/results/tx_quant/ledger.tsv
