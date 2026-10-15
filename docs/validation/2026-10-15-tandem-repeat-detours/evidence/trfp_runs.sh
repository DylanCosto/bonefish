#!/usr/bin/env bash
# Tandem-repeat false positives: sub-anchored gap fill on (dna preset) vs off (BONEFISH_EXP_OFF=subanchor)
# vs minimap2, on HG002 chr20 (PAW70337) and chr22 (PAW71238); BAMs kept for the detour count.
set -uo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
export BONEFISH="${BONEFISH:-$ROOT/bin/bonefish-trfp}"
R="$ROOT/results"
run() { # label chrom reads [env...]
  local label=$1 chrom=$2 reads=$3; shift 3
  [[ -s "$R/accuracy/$label/sv_accuracy.tsv" ]] && return
  echo "[$(date +%T)] $label"
  env "$@" KEEP_BAM=1 CHROM=$chrom READS="$reads" "$ROOT/bench/eval_accuracy.sh" "$label" "" > "$R/accuracy_$label.log" 2>&1
}
C20="$ROOT/data/reads/hg002_PAW70337_chr20.fq.gz"; C22="$ROOT/data/reads/hg002_PAW71238_chr22.fq.gz"
run trfp_off_chr20 chr20 "$C20" BONEFISH_EXP_OFF=subanchor
run trfp_on_chr22 chr22 "$C22"
run trfp_off_chr22 chr22 "$C22" BONEFISH_EXP_OFF=subanchor
run trfp_mm2_chr22 chr22 "$C22" ALIGNER=minimap2
echo "[$(date +%T)] done"
