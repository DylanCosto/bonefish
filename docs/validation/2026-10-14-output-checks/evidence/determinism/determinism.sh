#!/usr/bin/env bash
# Determinism and input-invariance: bonefish output must not depend on thread count, batch size,
# repeated runs, input encoding (FASTQ plain / gzip / bgzip / stdin / pipe / unaligned BAM / split
# into several files), index form (built on the fly from FASTA, .idx, raw) or SIMD code path.
# Each run's SAM is hashed without @PG (which records the command line). Also checks the
# built-in sorted BAM against the SAM output (same records) and the genome-wide HG001 sample
# against the published MD5 (v0.1.1 x86, reproduced on Apple M4 Pro).
# Usage: bench/determinism.sh OUTDIR     Env: BONEFISH
set -uo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
export PATH="${MAMBA_ROOT_PREFIX:-$HOME/micromamba}/envs/vareval/bin:$PATH"
B="${BONEFISH:-$ROOT/bonefish/target/release/bonefish}"
OUT="${1:?outdir}"; mkdir -p "$OUT"; W="$OUT/work"; mkdir -p "$W"
REF="$ROOT/data/ref/GRCh38_no_alt.fa"; RAW="$ROOT/data/index/bonefish.lrhq.raw"; IDX="$ROOT/data/index/bonefish.lrhq.idx"
Q="$ROOT/data/reads/hg001_wg_sub10.fq"
T="$OUT/results.tsv"; echo -e "group\trun\tmd5\twall_s" > "$T"

run() { # group, name, then the command (SAM on stdout)
  local g=$1 n=$2; shift 2
  local t0=$(date +%s.%N)
  local m; m=$("$@" 2> "$W/$n.log" | grep -v '^@PG' | md5sum | cut -d' ' -f1)
  printf "%s\t%s\t%s\t%.1f\n" "$g" "$n" "$m" "$(echo "$(date +%s.%N) - $t0" | bc)" | tee -a "$T"
}

# inputs derived from the same reads
[[ -s "$W/q.fq.gz" ]] || gzip -c "$Q" > "$W/q.fq.gz"
[[ -s "$W/q.bgz.fq.gz" ]] || bgzip -@8 -c "$Q" > "$W/q.bgz.fq.gz"
[[ -s "$W/q.u.bam" ]] || samtools import -@8 -0 "$Q" -O bam -o "$W/q.u.bam"
if [[ ! -s "$W/q.part2.fq" ]]; then
  n=$(( $(wc -l < "$Q") / 8 * 4 ))
  head -n "$n" "$Q" > "$W/q.part1.fq"; tail -n +$((n + 1)) "$Q" > "$W/q.part2.fq"
fi
cat "$RAW" > /dev/null

for t in 1 2 4 8 16; do run threads "t$t" "$B" -a -t "$t" "$RAW" "$Q"; done
for i in 1 2 3; do run repeat "t16_rep$i" "$B" -a -t 16 "$RAW" "$Q"; done
for k in 5M 50M 500M; do run batch "K$k" "$B" -a -t 16 -K "$k" "$RAW" "$Q"; done
run input fastq "$B" -a -t 16 "$RAW" "$Q"
run input gzip "$B" -a -t 16 "$RAW" "$W/q.fq.gz"
run input bgzip "$B" -a -t 16 "$RAW" "$W/q.bgz.fq.gz"
run input unaligned_bam "$B" -a -t 16 "$RAW" "$W/q.u.bam"
run input two_files "$B" -a -t 16 "$RAW" "$W/q.part1.fq" "$W/q.part2.fq"
run input stdin sh -c "cat '$Q' | '$B' -a -t 16 '$RAW' -"
run input stdin_bam sh -c "cat '$W/q.u.bam' | '$B' -a -t 16 '$RAW' -"
run input pipe bash -c "'$B' -a -t 16 '$RAW' <(zcat '$W/q.fq.gz')"
run index raw "$B" -a -t 16 "$RAW" "$Q"
run index idx "$B" -a -t 16 "$IDX" "$Q"
run index fasta "$B" -a -t 16 "$REF" "$Q"
for p in AVX2 SSE SCALAR; do run simd "force_$p" env "BONEFISH_FORCE_$p=1" "$B" -a -t 16 "$RAW" "$Q"; done
run simd auto "$B" -a -t 16 "$RAW" "$Q"

# sorted BAM: same records as the SAM output, and independent of threads
for t in 4 16; do
  "$B" -t "$t" -o "$W/sorted_t$t.bam" "$RAW" "$Q" 2> "$W/sorted_t$t.log"
  echo -e "sorted_bam\tt${t}_records\t$(samtools view "$W/sorted_t$t.bam" | md5sum | cut -d' ' -f1)\t" | tee -a "$T"
done
"$B" -a -t 16 "$RAW" "$Q" 2>/dev/null | grep -v '^@' | sort > "$W/sam.sorted.txt"
samtools view "$W/sorted_t16.bam" | sort > "$W/bam.sorted.txt"
# float tags (de:f) are compared as numbers: SAM text has 0.1750 where samtools prints 0.175
r=$(python3 - "$W/sam.sorted.txt" "$W/bam.sorted.txt" <<'PY'
import re, sys
def norm(path):
    return sorted("\t".join(t[:5] + repr(float(t[5:])) if re.match(r"^[A-Za-z][A-Za-z0-9]:f:", t) else t
                            for t in l.rstrip("\n").split("\t")) for l in open(path))
a, b = norm(sys.argv[1]), norm(sys.argv[2])
print(f"identical ({len(a)} records; f-type tags compared as numbers)" if a == b else f"DIFFERENT ({sum(1 for x, y in zip(a, b) if x != y)} records)")
PY
)
echo -e "sorted_bam\tsam_vs_bam_record_set\t$r\t" | tee -a "$T"
samtools view -H "$W/sorted_t16.bam" | grep '^@HD' | tee "$OUT/sorted_hd.txt"
rm -f "$W/sam.sorted.txt" "$W/bam.sorted.txt"

# genome-wide HG001 sample against the published MD5 (FASTA reference, as on the Mac)
run published hg001_wg_1x_t8_fasta "$B" -a -t 8 "$REF" "$ROOT/data/reads/hg001_wg_1x.fq.gz"
echo -e "published\texpected\tdeed7949d35f15bd8eed935be47a1d81\t" | tee -a "$T"

echo; echo "## distinct MD5s per group (1 = deterministic within the group)"
awk -F'\t' 'NR>1 && $1!="sorted_bam" && $1!="published" {print $1"\t"$3}' "$T" | sort -u | cut -f1 | uniq -c
