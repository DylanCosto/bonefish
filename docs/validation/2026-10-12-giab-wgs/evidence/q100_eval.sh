#!/usr/bin/env bash
# Re-score finished bench/sample_eval.sh HG002 Clair3 calls against GIAB HG002 v5.0q (T2T-Q100-based)
# small variants: overall, in the regions v5.0q adds over v4.2.1, and in GIAB v3.3 strata.
# Usage: bench/q100_eval.sh [OUT_DIR=data/wgs/HG002]       Writes OUT_DIR/{minimap2,bonefish}/q100.tsv
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
O="${1:-$ROOT/data/wgs/HG002}"; THREADS="${THREADS:-16}"
export PATH="${MAMBA_ROOT_PREFIX:-$HOME/micromamba}/envs/vareval/bin:$PATH"
T="$ROOT/data/truth"; ST="$T/strat"; SDF="$ROOT/data/ref/GRCh38_no_alt.sdf"
QV="$T/HG002_GRCh38_v5.0q_smvar.vcf.gz"; QB="$T/HG002_GRCh38_v5.0q_smvar.benchmark.bed"
OLD="$T/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed"
Q="$O/q100"; mkdir -p "$Q"
for kind in snps indels; do
  [[ -s "$Q/truth_$kind.vcf.gz.tbi" ]] || { bcftools view -v "$kind" -Oz -o "$Q/truth_$kind.vcf.gz" "$QV"; tabix -p vcf "$Q/truth_$kind.vcf.gz"; }
done
[[ -s "$Q/all.bed" ]] || sort -k1,1 -k2,2n "$QB" > "$Q/all.bed"
[[ -s "$Q/new_vs_v4.2.1.bed" ]] || "$ROOT/bench/bed_subtract.py" "$Q/all.bed" "$OLD" > "$Q/new_vs_v4.2.1.bed"
for s in alldifficultregions notinalldifficultregions segdups lowmappabilityall AllTandemRepeatsandHomopolymers_slop5; do
  [[ -s "$Q/$s.bed" ]] || "$ROOT/bench/bed_intersect.py" "$Q/all.bed" "$ST/GRCh38_$s.bed.gz" > "$Q/$s.bed"
done
for A in minimap2 bonefish; do
  D="$O/$A"; out="$D/q100.tsv"; echo -e "region\ttype\ttp\tfp\tfn\tprecision\trecall\tf1" > "$out.part"
  for r in all new_vs_v4.2.1 alldifficultregions notinalldifficultregions segdups lowmappabilityall AllTandemRepeatsandHomopolymers_slop5; do
    for kind in snps indels; do
      d="$Q/${A}_${r}_$kind"; rm -rf "$d"
      rtg vcfeval -b "$Q/truth_$kind.vcf.gz" -c "$D/calls_$kind.vcf.gz" -e "$Q/$r.bed" -t "$SDF" -o "$d" --threads "$THREADS" > /dev/null
      awk -v r="$r" -v k="$kind" '$1=="None"{printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", r, k, $2, $4, $5, $6, $7, $8}' "$d/summary.txt" >> "$out.part"
      [[ $r == all ]] || rm -rf "$d"
    done
  done
  mv "$out.part" "$out"
done
paste <(cut -f1-2,8 "$O/minimap2/q100.tsv") <(cut -f3-5,8 "$O/minimap2/q100.tsv") <(cut -f3-5,8 "$O/bonefish/q100.tsv") | column -t
