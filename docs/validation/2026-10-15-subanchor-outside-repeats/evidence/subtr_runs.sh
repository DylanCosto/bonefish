#!/usr/bin/env bash
# Validation of the tandem-repeat guard on sub-anchoring (dna preset, --gap-subanchor-tr 150):
# chr20/chr22 accuracy (SNP/indel/SV), GIAB tandem-repeat benchmark + detours, simulated reads.
set -uo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"; R="$ROOT/results"
export BONEFISH="$ROOT/bin/bonefish-subtr"
for c in chr20 chr22; do
  label=subtr150_$c; [[ -s "$R/accuracy/$label/sv_accuracy.tsv" ]] && continue
  reads="$ROOT/data/reads/hg002_PAW70337_chr20.fq.gz"; [[ $c == chr22 ]] && reads="$ROOT/data/reads/hg002_PAW71238_chr22.fq.gz"
  echo "[$(date +%T)] $label"
  KEEP_BAM=1 CHROM=$c READS="$reads" "$ROOT/bench/eval_accuracy.sh" "$label" "" > "$R/accuracy_$label.log" 2>&1
done
echo "[$(date +%T)] tandem-repeat scoring"
"$ROOT/bench/trfp_score.sh" chr20 subtr150_chr20 > "$R/subtr/score20.log" 2>&1
"$ROOT/bench/trfp_score.sh" chr22 subtr150_chr22 > "$R/subtr/score22.log" 2>&1
echo "[$(date +%T)] simulated reads"
mkdir -p "$R/sim_wg_subtr" && cp "$R"/sim_wg/*minimap2* "$R/sim_wg_subtr/"
BONEFISH="$BONEFISH" "$ROOT/bench/sim_wg_eval.sh" "$R/sim_wg_subtr" > "$R/sim_wg_subtr.log" 2>&1
echo "[$(date +%T)] done"
