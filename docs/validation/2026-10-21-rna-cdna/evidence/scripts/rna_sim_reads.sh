#!/usr/bin/env bash
# Simulated spliced reads with known intron chains (as run on 21 Oct 2026): 30,000 random
# multi-exon GENCODE v49 transcripts, Badread 0.4 reads at cDNA-like (98% mean identity, both
# strands) and direct-RNA-like (95%, sense strand only) accuracy, no junk/chimeras/adapters.
set -euo pipefail
R=${BENCH}; D=$R/data/rna/sim; mkdir -p $D; cd $D
B=$HOME/micromamba/envs/sim/bin/badread
python3 $R/bench/rna_sim_transcripts.py $R/data/ref/GRCh38_no_alt.fa ../gencode.v49.bed12 30000 7 tx
common=(--reference tx.fa --length 10000,5000 --junk_reads 0 --random_reads 0 --chimeras 0
        --glitches 0,0,0 --start_adapter_seq "" --end_adapter_seq "")
$B simulate "${common[@]}" --quantity 80M --identity 98,100,1.5 --seed 11 2>cdna.log | gzip -1 > sim_cdna.fq.gz
$B simulate "${common[@]}" --quantity 160M --identity 95,99,3 --seed 12 2>drna.log \
  | python3 -c "
import sys
lines = sys.stdin.readlines()
for i in range(0, len(lines), 4):
    if ',+strand,' in lines[i]: sys.stdout.writelines(lines[i:i + 4])" | gzip -1 > sim_drna.fq.gz
