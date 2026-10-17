#!/usr/bin/env bash
# SIMD parity and speed check for bonefish on one x86_64 machine (written for an AVX-512 test).
#
# Maps the genome-wide HG001 sample with every SIMD path the CPU supports and checks that each
# gives the same output, and the same output as the reference x86 (AVX2) and arm64 runs:
#   v0.1.5-v0.1.6: MD5 9bc7ff3e9f498cc1530b2b1d78c47e10   (v0.1.1-v0.1.3: deed7949d35f15bd8eed935be47a1d81)
# MD5 is of the SAM with @PG lines removed. Then times the default (AVX-512 where present) against
# BONEFISH_FORCE_AVX2=1, alternating, and prints a report to paste back.
#
# Usage: scripts/validation/simd_parity.sh BONEFISH REFERENCE.fa[.gz] hg001_wg_1x.fq.gz [WORKDIR]
# Env: THREADS (all cores), REPEATS (timing rounds, 3)
set -euo pipefail
B=$(readlink -f "$1"); REF=$(readlink -f "$2"); READS=$(readlink -f "$3"); W=${4:-simd_parity_out}
THREADS=${THREADS:-$(nproc)}; REPEATS=${REPEATS:-3}
EXPECT=9bc7ff3e9f498cc1530b2b1d78c47e10
mkdir -p "$W"; cd "$W"
R=report.txt; : > $R
say() { echo "$*" | tee -a $R; }

say "== machine"
say "$(lscpu | grep -E '^(Model name|CPU\(s\):|Thread|Socket|NUMA node\(s\)|Hypervisor vendor)' | sed 's/  */ /g')"
say "memory: $(free -g | awk '/Mem:/{print $2}') GB; threads used: $THREADS"
flags=$(grep -m1 '^flags' /proc/cpuinfo)
for f in avx2 avx512f avx512bw avx512vl avx512_vnni avx512_bf16; do
  [[ $flags == *" $f"* ]] && say "  $f: yes" || say "  $f: no"
done
[[ $flags == *" avx512bw"* ]] || say "WARNING: the guest does not see avx512bw, so the AVX-512 paths cannot run (check the VM's CPU model / passthrough)."
say "bonefish: $("$B" --version)  sha256 $(sha256sum "$B" | cut -c1-16)"
say "reads sha256: $(sha256sum "$READS" | cut -c1-16) (expected 2edff9a3834d5d55)"

say "== index"
IDX=grch38.raw
if [[ ! -s $IDX ]]; then
  /usr/bin/time -f "%e s wall, %M KB" -o index.time "$B" -t "$THREADS" -d $IDX --idx-raw "$REF" 2> index.log
  say "built in $(cat index.time)"
fi
cat $IDX > /dev/null   # page it in so the first timed run is not slower

# name|environment
PATHS=("default (AVX-512 where present)|" "AVX2 (BONEFISH_FORCE_AVX2=1)|BONEFISH_FORCE_AVX2=1" "SSE4.1 (BONEFISH_FORCE_SSE=1)|BONEFISH_FORCE_SSE=1")
run() { # env -> prints "wall cpu md5"
  env $1 /usr/bin/time -f "%e %U %S" -o t.tmp "$B" -a -t "$THREADS" $IDX "$READS" 2>/dev/null | grep -v '^@PG' | md5sum | cut -d' ' -f1 > m.tmp
  read -r w u s < t.tmp; echo "$w $(echo "$u + $s" | bc) $(cat m.tmp)"
}

say "== output identity (expected $EXPECT)"
ok=1
for p in "${PATHS[@]}"; do
  name=${p%%|*}; e=${p#*|}
  read -r w c m <<< "$(run "$e")"
  [[ $m == "$EXPECT" ]] && st=MATCH || { st=DIFFERS; ok=0; }
  say "  $name: $m $st (${w} s wall, ${c} CPU-s)"
done
[[ $ok == 1 ]] && say "RESULT: every SIMD path reproduces the reference output" || say "RESULT: OUTPUT DIFFERS - send report.txt"

say "== speed: default vs AVX2, $REPEATS alternating rounds"
for i in $(seq 1 "$REPEATS"); do
  for p in "${PATHS[@]:0:2}"; do
    name=${p%%|*}; e=${p#*|}
    read -r w c m <<< "$(run "$e")"
    say "  round $i $name: ${w} s wall, ${c} CPU-s"
  done
done
rm -f t.tmp m.tmp
say "== done; report in $(pwd)/$R"
