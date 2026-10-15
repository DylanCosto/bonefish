#!/usr/bin/env bash
# Mapping accuracy on the whole-genome Badread sets (bench/sim_wg.sh): bonefish and minimap2 2.31
# (lr:hq; also map-ont for the low-accuracy set) against GRCh38 no-alt, scored by bench/sim_eval.py.
# Usage: bench/sim_wg_eval.sh OUTDIR     Env: BONEFISH
set -uo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
B="${BONEFISH:-$ROOT/bonefish/target/release/bonefish}"
MM="${MAMBA_ROOT_PREFIX:-$HOME/micromamba}/envs/ontaln/bin/minimap2"
RAW="$ROOT/data/index/bonefish.lrhq.raw"; MMI="$ROOT/data/index/mm2.lrhq.mmi"; REF="$ROOT/data/ref/GRCh38_no_alt.fa"
S="$ROOT/data/sim/wg"; OUT="${1:?outdir}"; mkdir -p "$OUT"
for p in r10 lowacc short ul; do
  q="$S/$p.fq.gz"
  runs=("bonefish|$B -a -t 16 $RAW" "minimap2_lrhq|$MM -ax lr:hq -t 16 $MMI")
  [[ $p == lowacc ]] && runs+=("minimap2_mapont|$MM -ax map-ont -t 16 $REF")
  for r in "${runs[@]}"; do
    name=${r%%|*}; cmd=${r#*|}
    [[ -s "$OUT/$p.$name.md" ]] && continue
    /usr/bin/time -f "%e s wall, %U+%S s CPU, %M KB peak" -o "$OUT/$p.$name.time" $cmd "$q" 2> "$OUT/$p.$name.log" \
      | python3 "$ROOT/bench/sim_eval.py" "$p $name" "$q" "$OUT/$p.$name.curve.tsv" > "$OUT/$p.$name.md"
    echo "time: $(cat "$OUT/$p.$name.time")" >> "$OUT/$p.$name.md"
    cat "$OUT/$p.$name.md"
  done
done
