#!/usr/bin/env bash
# Thread scaling and memory on the genome-wide HG001 sample (161,762 reads, 2.7 Gbp; SAM output
# hashed, no disk I/O): bonefish at 1-16 threads, minimap2 2.31 lr:hq at 4/8/16 threads; then
# peak memory for small-machine settings (fewer threads, smaller batches, split index).
# Run on an otherwise idle machine. Usage: bench/scaling.sh OUTDIR     Env: BONEFISH
set -uo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
B="${BONEFISH:-$ROOT/bonefish/target/release/bonefish}"
MM="${MAMBA_ROOT_PREFIX:-$HOME/micromamba}/envs/ontaln/bin/minimap2"
RAW="$ROOT/data/index/bonefish.lrhq.raw"; IDX="$ROOT/data/index/bonefish.lrhq.idx"; MMI="$ROOT/data/index/mm2.lrhq.mmi"
REF="$ROOT/data/ref/GRCh38_no_alt.fa"; Q="$ROOT/data/reads/hg001_wg_1x.fq.gz"
OUT="${1:?outdir}"; mkdir -p "$OUT"; T="$OUT/scaling.tsv"
[[ -s "$T" ]] || echo -e "aligner\tsetting\tthreads\twall_s\tcpu_s\tpeak_rss_gb\tmd5" > "$T"
run() { # aligner, setting, threads, then command
  local a=$1 s=$2 t=$3; shift 3
  grep -q "^$a	$s	$t	" "$T" && return
  /usr/bin/time -f "%e %U %S %M" -o "$OUT/time.tmp" "$@" "$Q" 2> "$OUT/$a.$s.t$t.log" | grep -v '^@PG' | md5sum | cut -c1-12 > "$OUT/md5.tmp"
  read -r w u sy m < "$OUT/time.tmp"
  printf "%s\t%s\t%s\t%.1f\t%.0f\t%.2f\t%s\n" "$a" "$s" "$t" "$w" "$(echo "$u + $sy" | bc)" "$(echo "scale=2; $m / 1048576" | bc)" "$(cat "$OUT/md5.tmp")" | tee -a "$T"
}
cat "$RAW" > /dev/null
for t in 16 12 8 4 2 1; do run bonefish raw_index "$t" "$B" -a -t "$t" "$RAW"; done
cat "$MMI" > /dev/null
for t in 16 8 4; do run minimap2 lrhq "$t" "$MM" -ax lr:hq -t "$t" "$MMI"; done
# small-machine settings
cat "$IDX" > /dev/null
run bonefish idx_index 4 "$B" -a -t 4 "$IDX"
run bonefish raw_index_K10M 4 "$B" -a -t 4 -K 10M "$RAW"
run bonefish idx_index_K10M 4 "$B" -a -t 4 -K 10M "$IDX"
run bonefish split_index_I2G 4 "$B" -a -t 4 -K 10M -I 2G --split-prefix "$OUT/split.tmp" "$REF"
run bonefish split_index_I1G 4 "$B" -a -t 4 -K 10M -I 1G --split-prefix "$OUT/split.tmp" "$REF"
rm -f "$OUT"/time.tmp "$OUT"/md5.tmp "$OUT"/split.tmp*
