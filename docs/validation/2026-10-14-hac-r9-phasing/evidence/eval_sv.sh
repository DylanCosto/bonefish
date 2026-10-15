#!/usr/bin/env bash
# Structural-variant check for a Phase 3 evaluation: Sniffles2 on results/accuracy/LABEL/aln.bam,
# Truvari against GIAB HG002 v5.0q GRCh38 SVs (>= 50 bp) within the benchmark regions of CHROM
# (default chr20).
# Matching ignores sequence similarity (--pctseq 0): size and position only, since Sniffles
# deletion/insertion sequences are not always populated.
#
# Usage: bench/eval_sv.sh LABEL
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
LABEL="${1:?label}"
CHROM="${CHROM:-chr20}"
THREADS="${THREADS:-16}"
export PATH="${MAMBA_ROOT_PREFIX:-$HOME/micromamba}/envs/vareval/bin:$PATH"
REF="$ROOT/data/ref/GRCh38_no_alt.fa"
T="$ROOT/data/truth"
OUT="$ROOT/results/accuracy/$LABEL"
BAM="$OUT/aln.bam"
[[ -s "$BAM" ]] || { echo "no BAM at $BAM"; exit 1; }

TV="$T/HG002_${CHROM}_v5.0q_stvar.vcf.gz"; TB="$T/HG002_${CHROM}_v5.0q_stvar.bed"
if [[ ! -s "$TV.tbi" ]]; then
  bcftools view -r "$CHROM" -Oz -o "$TV" "$T/HG002_GRCh38_v5.0q_stvar.vcf.gz"
  tabix -p vcf "$TV"
  awk -v c="$CHROM" '$1==c' "$T/HG002_GRCh38_v5.0q_stvar.benchmark.bed" > "$TB"
fi

SV="$OUT/sniffles.vcf.gz"
if [[ ! -s "$SV" ]]; then
  sniffles --input "$BAM" --vcf "$OUT/sniffles.vcf" --reference "$REF" --threads "$THREADS" \
    --regions "$TB" > "$OUT/sniffles.log" 2>&1
  bcftools sort -Oz -o "$SV" "$OUT/sniffles.vcf" && tabix -p vcf "$SV" && rm -f "$OUT/sniffles.vcf"
fi

rm -rf "$OUT/truvari"
truvari bench -b "$TV" -c "$SV" --includebed "$TB" --passonly --pctseq 0 -o "$OUT/truvari" > "$OUT/truvari.log" 2>&1
python3 - "$LABEL" "$OUT/truvari/summary.json" <<'EOF' | tee "$OUT/sv_accuracy.tsv"
import json, sys
s = json.load(open(sys.argv[2]))
print("label\tTP-base\tFP\tFN\tprecision\trecall\tf1")
print(f"{sys.argv[1]}\t{s['TP-base']}\t{s['FP']}\t{s['FN']}\t{s['precision']:.4f}\t{s['recall']:.4f}\t{s['f1']:.4f}")
EOF
