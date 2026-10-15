#!/usr/bin/env bash
# Detour repair variant (BONEFISH_EXP_SUB_REPAIR=30) on chr20 and chr22, after trfp_runs.sh.
set -uo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
R="$ROOT/results"
while ps aux | grep -q '[t]rfp_runs\.sh'; do sleep 20; done
export BONEFISH="$ROOT/bin/bonefish-trfp2"
for c in chr20 chr22; do
  label=trfp_rep30_$c
  [[ -s "$R/accuracy/$label/sv_accuracy.tsv" ]] && continue
  reads="$ROOT/data/reads/hg002_PAW70337_chr20.fq.gz"; [[ $c == chr22 ]] && reads="$ROOT/data/reads/hg002_PAW71238_chr22.fq.gz"
  echo "[$(date +%T)] $label"
  BONEFISH_EXP_SUB_REPAIR=30 KEEP_BAM=1 CHROM=$c READS="$reads" "$ROOT/bench/eval_accuracy.sh" "$label" "" > "$R/accuracy_$label.log" 2>&1
done
echo "[$(date +%T)] done"
