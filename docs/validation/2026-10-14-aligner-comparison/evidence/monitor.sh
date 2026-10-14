#!/usr/bin/env bash
# Resource sampler for benchmark runs.
# Appends one CSV row every INTERVAL seconds to results/logs/resources.csv and prints:
#   - a one-line summary every SUMMARY_EVERY samples
#   - an ALERT line immediately on low memory, swap growth, low disk, or idle CPU
#     while an aligner is running (a sign of an I/O or threading bottleneck)
#
# Usage: bench/monitor.sh [interval_s] [summary_every]
# Env: MON_CSV (output CSV, default results/logs/resources.csv), MON_STAGE (a file whose first line
#      names the current pipeline stage; recorded in a "stage" column)
set -uo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
INTERVAL="${1:-5}"
SUMMARY_EVERY="${2:-36}"
CSV="${MON_CSV:-$ROOT/results/logs/resources.csv}"
STAGE_FILE="${MON_STAGE:-}"
PROCS='minimap2|bonefish|rammap|dorado|mm2plus|winnowmap|ngmlr|lra|pbrun|pbminimap2|samtools|perf|python|pypy|longphase|whatshap|sniffles|java|aws|bgzip|parallel'
mkdir -p "$(dirname "$CSV")"
[[ -s "$CSV" ]] || echo "time,cpu_pct,mem_used_gb,mem_avail_gb,swap_used_gb,disk_free_gb,proc,proc_cpu_pct,proc_rss_gb,gpu_util_pct,gpu_mem_gb,stage,iowait_pct" > "$CSV"

read_cpu() { awk '/^cpu /{print $2+$3+$4+$5+$6+$7+$8, $5+$6, $6}' /proc/stat; }
read -r t0 i0 w0 < <(read_cpu)
swap0=""
n=0 sum_cpu=0 peak_rss=0 peak_mem=0 min_avail=999 low_cpu=0

while true; do
  sleep "$INTERVAL"
  read -r t1 i1 w1 < <(read_cpu)
  cpu=$(( (100 * ((t1 - t0) - (i1 - i0))) / ((t1 - t0) > 0 ? (t1 - t0) : 1) ))
  iow=$(( (100 * (w1 - w0)) / ((t1 - t0) > 0 ? (t1 - t0) : 1) ))
  t0=$t1 i0=$i1 w0=$w1
  stage=$( [[ -n "$STAGE_FILE" && -s "$STAGE_FILE" ]] && head -1 "$STAGE_FILE" || echo -)

  read -r used avail swap < <(free -m | awk '/^Mem:/{u=$3; a=$7} /^Swap:/{s=$3} END{printf "%.1f %.1f %.2f\n", u/1024, a/1024, s/1024}')
  disk=$(df -BG --output=avail "$ROOT/data/" | tail -1 | tr -dc '0-9')   # the data volume
  # Largest benchmark process by RSS, with its CPU share across all its threads.
  read -r pname pcpu prss < <(ps -eo comm=,pcpu=,rss= | awk -v re="^($PROCS)" '$1 ~ re {if ($3 > r) {r=$3; c=$2; n=$1}} END{if (n) printf "%s %s %.1f\n", n, c, r/1048576; else print "- 0 0"}')
  read -r gutil gmem < <(nvidia-smi --query-gpu=utilization.gpu,memory.used --format=csv,noheader,nounits 2>/dev/null | awk -F', ' '{printf "%s %.1f\n", $1, $2/1024}')

  echo "$(date '+%F %T'),$cpu,$used,$avail,$swap,$disk,$pname,$pcpu,$prss,${gutil:-0},${gmem:-0},$stage,$iow" >> "$CSV"

  n=$((n + 1)); sum_cpu=$((sum_cpu + cpu))
  awk -v a="$prss" -v b="$peak_rss" 'BEGIN{exit !(a > b)}' && peak_rss=$prss
  awk -v a="$used" -v b="$peak_mem" 'BEGIN{exit !(a > b)}' && peak_mem=$used
  awk -v a="$avail" -v b="$min_avail" 'BEGIN{exit !(a < b)}' && min_avail=$avail
  swap0=${swap0:-$swap}

  awk -v a="$avail" 'BEGIN{exit !(a < 2.0)}' && echo "ALERT $(date +%T) low memory: ${avail} GB available ($pname RSS ${prss} GB)"
  awk -v s="$swap" -v s0="$swap0" 'BEGIN{exit !(s - s0 > 0.5)}' && { echo "ALERT $(date +%T) swap grew to ${swap} GB"; swap0=$swap; }
  (( disk < 15 )) && echo "ALERT $(date +%T) low disk: ${disk} GB free"
  # Idle CPU while an aligner runs: alert after 6 consecutive samples (~30 s). Index loading is
  # single-threaded for ~15-25 s and should not trip it.
  if [[ "$pname" =~ ^(minimap2|bonefish|rammap|pbminimap2)$ ]] && (( cpu < 50 )); then
    low_cpu=$((low_cpu + 1)); (( low_cpu == 6 )) && echo "ALERT $(date +%T) CPU only ${cpu}% while $pname runs (possible I/O or serial bottleneck)"
  else
    low_cpu=0
  fi

  if (( n % SUMMARY_EVERY == 0 )); then
    echo "$(date +%T) avg CPU $((sum_cpu / SUMMARY_EVERY))% | mem peak ${peak_mem} GB, min avail ${min_avail} GB | swap ${swap} GB | disk ${disk} GB | top: $pname ${pcpu}% CPU, peak RSS ${peak_rss} GB | GPU ${gutil:-0}% ${gmem:-0} GB"
    sum_cpu=0 peak_rss=0 peak_mem=0 min_avail=999
  fi
done
