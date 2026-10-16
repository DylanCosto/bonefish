#!/usr/bin/env bash
# Whole-genome re-validation of a new bonefish build against the kept minimap2 results of the
# 12 Oct run (data/wgs/SAMPLE/minimap2, linked): bench/sample_eval.sh with SUFFIX=_v2, then for
# HG002 the GIAB v5.0q small-variant scoring and the CMRG / tandem-repeat benchmarks.
# Reads come from bench/wgs_fetch.sh (data/wgs/SAMPLE_v2/reads.fq.gz). Usage: BONEFISH=... bench/wgs_v2.sh
set -uo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"; W="$ROOT/data/wgs"
: "${BONEFISH:?set BONEFISH}"; export BONEFISH
# Samples are processed in the order their reads become ready.
todo=(HG002 HG005)
while (( ${#todo[@]} )); do
  S=""
  for c in "${todo[@]}"; do [[ -s "$W/${c}_v2/reads.fq.gz" ]] && { S=$c; break; }; done
  if [[ -z $S ]]; then
    # drop samples with neither reads nor a fetch running or queued
    keep=(); for c in "${todo[@]}"; do
      if pgrep -f "wgs_fetch.sh $c" > /dev/null || [[ -s "$W/${c}_v2/reads.fq.gz" ]]; then keep+=("$c")
      else echo "[$(date '+%F %T')] $c: no reads and no fetch running, skipped"; fi
    done; todo=("${keep[@]}"); sleep 60; continue
  fi
  todo=($(printf '%s\n' "${todo[@]}" | grep -vx "$S"))
  O="$W/${S}_v2"
  [[ -e "$O/minimap2" ]] || ln -s "$W/$S/minimap2" "$O/minimap2"
  for f in "$W/$S"/truth_* "$W/$S"/bench_*.bed; do [[ -e "$O/$(basename "$f")" ]] || cp "$f" "$O/"; done
  echo "[$(date '+%F %T')] $S: sample_eval"
  SUFFIX=_v2 "$ROOT/bench/sample_eval.sh" "$S" > "$O/sample_eval.out" 2>&1 || echo "[$(date '+%F %T')] $S: sample_eval failed"
  if [[ $S == HG002 ]]; then
    echo "[$(date '+%F %T')] HG002: v5.0q"
    [[ -e "$O/q100" ]] || cp -r "$W/HG002/q100" "$O/q100" 2>/dev/null; rm -rf "$O"/q100/bonefish_* 
    SKIP_DONE=1 "$ROOT/bench/q100_eval.sh" "$O" > "$O/q100.out" 2>&1
    echo "[$(date '+%F %T')] HG002: CMRG + tandem repeats"
    WGS_DIR="$O" ALIGNERS=bonefish "$ROOT/bench/cmrg_tr_eval.sh" "$ROOT/results/cmrg_tr_v2" > "$ROOT/results/cmrg_tr_v2.log" 2>&1
  fi
done
echo "[$(date '+%F %T')] done"
