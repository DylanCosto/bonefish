#!/usr/bin/env bash
# Score results/accuracy runs on the GIAB tandem-repeat benchmark (Tier 1 regions of one chromosome;
# Clair3 PASS + Sniffles2 PASS merged, multi-allelic split; truvari bench --sizemin 5 --pick ac, then
# refine) and count deletion-insertion detours in Tier 1 repeats (bench/detours.py).
# Usage: bench/trfp_score.sh CHROM LABEL...     Output: results/trfp/score_<CHROM>.tsv
set -uo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
export PATH="${MAMBA_ROOT_PREFIX:-$HOME/micromamba}/envs/vareval/bin:$PATH"
C=$1; shift; O="$ROOT/results/trfp"; TR="$ROOT/data/truth/tr"; REF="$ROOT/data/ref/GRCh38_no_alt.fa"
[[ -s "$O/tier1_$C.bed" ]] || zcat "$TR/HG002_GRCh38_TandemRepeats_v1.0.1_Tier1.bed.gz" | awk -v c="$C" '$1==c' > "$O/tier1_$C.bed"
T="$O/score_$C.tsv"; [[ -s "$T" ]] || echo -e "label\tbench_tp\tbench_fp\tbench_fn\tbench_f1\tfp_ge50\tfn_ge50\trefine_tp\trefine_fp\trefine_fn\trefine_f1\tdetour_reads\tdetour_loci" > "$T"
for L in "$@"; do
  A="$ROOT/results/accuracy/$L"; W="$O/$L"; mkdir -p "$W"
  echo SAMPLE > "$W/s.txt"
  bcftools view -f PASS "$A/clair3/merge_output.vcf.gz" | bcftools reheader -s "$W/s.txt" | bcftools view -Oz -o "$W/c.vcf.gz"
  bcftools view -f PASS "$A/sniffles.vcf.gz" | bcftools reheader -s "$W/s.txt" | bcftools view -Oz -o "$W/v.vcf.gz"
  tabix -f "$W/c.vcf.gz"; tabix -f "$W/v.vcf.gz"
  bcftools concat -a "$W/c.vcf.gz" "$W/v.vcf.gz" 2>/dev/null | bcftools norm -m -any 2>/dev/null | bcftools sort -Oz -o "$W/merged.vcf.gz" 2>/dev/null
  tabix -f "$W/merged.vcf.gz"; rm -f "$W"/c.vcf.gz* "$W"/v.vcf.gz* "$W/s.txt"
  rm -rf "$W/bench"
  truvari bench -b "$TR/HG002_GRCh38_TandemRepeats_v1.0.1.vcf.gz" -c "$W/merged.vcf.gz" --includebed "$O/tier1_$C.bed" \
    --sizemin 5 --pick ac -o "$W/bench" > "$W/bench.log" 2>&1
  truvari refine --use-original-vcfs --reference "$REF" --threads 8 "$W/bench" > "$W/refine.log" 2>&1
  fp50=$(bcftools view -H -i 'abs(strlen(ALT)-strlen(REF))>=50' "$W/bench/fp.vcf.gz" | wc -l)
  fn50=$(bcftools view -H -i 'abs(strlen(ALT)-strlen(REF))>=50' "$W/bench/fn.vcf.gz" | wc -l)
  det=$(python3 "$ROOT/bench/detours.py" "$A/aln.bam" "$O/tier1_$C.bed" "$C" | awk -F'\t' '{split($4,a," "); split($5,b," "); print a[3]"\t"b[2]}')
  python3 - "$W/bench/summary.json" "$W/bench/refine.variant_summary.json" "$L" "$fp50" "$fn50" "$det" <<'PY' >> "$T"
import json, sys
b = json.load(open(sys.argv[1]))
try: r = json.load(open(sys.argv[2]))
except Exception: r = {"TP-base": "", "FP": "", "FN": "", "f1": float("nan")}
print("\t".join(map(str, [sys.argv[3], b["TP-base"], b["FP"], b["FN"], round(b["f1"], 4), sys.argv[4], sys.argv[5],
                          r["TP-base"], r["FP"], r["FN"], round(r["f1"], 4), sys.argv[6]])))
PY
done
column -t -s$'\t' "$T"
