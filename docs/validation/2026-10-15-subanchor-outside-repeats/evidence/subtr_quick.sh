#!/usr/bin/env bash
# Quick screen for the tandem-repeat guard on sub-anchoring: align HG002 chr22 (27x) with one config and
# count deletion-insertion detours in Tier 1 repeats. Usage: bench/subtr_quick.sh LABEL "ENV=.." "ARGS"
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"; L=$1; ENVS=$2; ARGS=$3
B=${BONEFISH:-$ROOT/bonefish/target/release/bonefish}; O=$ROOT/results/subtr/$L; mkdir -p "$O"
export PATH=$HOME/micromamba/envs/vareval/bin:$PATH
env $ENVS /usr/bin/time -f "%e %U %S %M" -o "$O/align.time" "$B" -t 16 $ARGS -o "$O/aln.bam" \
  "$ROOT/data/index/bonefish.lrhq.raw" "$ROOT/data/reads/${READS_CHR:-hg002_PAW71238_chr22}.fq.gz" 2> "$O/align.log"
det=$(python3 "$ROOT/bench/detours.py" "$O/aln.bam" "$ROOT/results/trfp/tier1_${CHR:-chr22}.bed" ${CHR:-chr22})
read w u s m < "$O/align.time"
printf "%s\t%s CPU-s\t%s\n" "$L" "$(echo "$u+$s" | bc)" "$det" | tee -a "$ROOT/results/subtr/quick.tsv"
