#!/usr/bin/env bash
# Whole-genome Badread simulations for mapping-accuracy checks (bench/sim_eval.py).
# Four read profiles, each simulated per chromosome (chr1-22, X, Y of GRCh38 no-alt) at the same
# coverage so reads are spread over the genome in proportion to length:
#   r10    R10.4.1 sup-like: identity 98/99.5/1.5 (mean/max/sd), length 15/13 kb, 0.33x
#   lowacc HAC / R9-like:    identity 93/98/4,     length 15/13 kb, 0.16x
#   short  cfDNA/amplicon:   identity 98/99.5/1.5, length 0.5/0.3 kb, 0.033x
#   ul     ultra-long:       identity 98/99.5/1.5, length 80/50 kb, 0.2x
# Badread defaults otherwise (nanopore2023 error/qscore models, ~1% each of junk, random and chimeric
# reads, adapters, glitches). Output: data/sim/wg/<profile>.fq.gz. Env: JOBS (parallel Badread runs).
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
S="$ROOT/data/sim/wg"; BR="${MAMBA_ROOT_PREFIX:-$HOME/micromamba}/envs/sim/bin/badread"
export PATH="${MAMBA_ROOT_PREFIX:-$HOME/micromamba}/envs/vareval/bin:$PATH"
JOBS="${JOBS:-10}"
mkdir -p "$S/ref" "$S/parts"
CHROMS=$(awk '$1 ~ /^chr([0-9]+|X|Y)$/ {print $2, $1}' "$ROOT/data/ref/GRCh38_no_alt.fa.fai" | sort -k1,1nr | cut -d' ' -f2)
for c in $CHROMS; do [[ -s "$S/ref/$c.fa" ]] || samtools faidx "$ROOT/data/ref/GRCh38_no_alt.fa" "$c" > "$S/ref/$c.fa"; done

declare -A ARGS=(
  [r10]="--identity 98,99.5,1.5 --length 15000,13000 --quantity 0.33x"
  [lowacc]="--identity 93,98,4 --length 15000,13000 --quantity 0.16x"
  [short]="--identity 98,99.5,1.5 --length 500,300 --quantity 0.033x"
  [ul]="--identity 98,99.5,1.5 --length 80000,50000 --quantity 0.2x"
)
jobs=()
seed=0
for p in r10 lowacc short ul; do
  for c in $CHROMS; do
    seed=$((seed + 1))
    [[ -s "$S/parts/$p.$c.fq.gz" ]] || jobs+=("$p $c $seed")
  done
done
# bash can't export associative arrays: pass each profile's arguments through the environment
for p in r10 lowacc short ul; do export "ARGS_$p=${ARGS[$p]}"; done
run_one() {
  set -euo pipefail
  read -r p c seed <<< "$1"
  local a="ARGS_$p"
  "$BR" simulate --reference "$S/ref/$c.fa" --error_model nanopore2023 --qscore_model nanopore2023 \
    ${!a} --seed "$seed" 2> /dev/null | gzip -1 > "$S/parts/$p.$c.fq.gz.tmp"
  mv "$S/parts/$p.$c.fq.gz.tmp" "$S/parts/$p.$c.fq.gz"
  echo "[$(date +%T)] $p $c done"
}
export -f run_one; export S BR
printf '%s\n' "${jobs[@]}" | xargs -P "$JOBS" -I{} bash -c 'run_one "$@"' _ {}
for p in r10 lowacc short ul; do
  [[ -s "$S/$p.fq.gz" ]] || { for c in $CHROMS; do cat "$S/parts/$p.$c.fq.gz"; done > "$S/$p.fq.gz"; }
  zcat "$S/$p.fq.gz" | awk -v p="$p" 'NR%4==2{n++;b+=length($0)}END{printf "%s: %d reads, %.1f Mbp\n", p, n, b/1e6}'
done
