#!/usr/bin/env bash
# Phasing and haplotagging from each aligner's chr20 BAM (HG002 PAW70337, 55x, kept by
# eval_accuracy.sh with KEEP_BAM=1): longphase (as in ONT's wf-human-variation) phases the Clair3
# PASS SNPs and indels; whatshap compare scores the phased SNVs against GIAB v5.0q (Q100
# assembly, fully phased), restricted to its benchmark regions; whatshap stats gives block N50;
# longphase haplotag gives the fraction of reads assigned to a haplotype.
# Usage: bench/phasing.sh OUTDIR LABEL...      (LABELs of results/accuracy runs)
set -uo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
export PATH="${MAMBA_ROOT_PREFIX:-$HOME/micromamba}/envs/vareval/bin:$PATH"
OUT="${1:?outdir}"; shift; mkdir -p "$OUT"
REF="$ROOT/data/ref/GRCh38_no_alt.fa"; T="$ROOT/data/truth"
TV="$OUT/truth_chr20.vcf.gz"
if [[ ! -s "$TV.tbi" ]]; then
  awk '$1=="chr20"' "$T/HG002_GRCh38_v5.0q_smvar.benchmark.bed" > "$OUT/truth_chr20.bed"
  bcftools view -r chr20 -T "$OUT/truth_chr20.bed" -i 'FILTER=="." || FILTER=="PASS"' -Oz -o "$TV" "$T/HG002_GRCh38_v5.0q_smvar.vcf.gz"
  tabix -p vcf "$TV"
fi
echo -e "label\tphased_snvs_compared\tswitch_errors\tswitch_error_rate\thamming_rate\tblock_n50_bp\tblocks\treads_haplotagged" > "$OUT/summary.tsv"
for label in "$@"; do
  A="$ROOT/results/accuracy/$label"; O="$OUT/$label"; mkdir -p "$O"
  bcftools view -f PASS -Oz -o "$O/calls.vcf.gz" "$A/clair3/merge_output.vcf.gz"; tabix -f -p vcf "$O/calls.vcf.gz"
  longphase phase -s "$O/calls.vcf.gz" -b "$A/aln.bam" -r "$REF" -t 16 --ont --indels -o "$O/phased" > "$O/longphase.log" 2>&1
  bgzip -f "$O/phased.vcf"; tabix -f -p vcf "$O/phased.vcf.gz"
  # whatshap compares a sample present in both files: give the calls the truth's sample name
  bcftools query -l "$TV" > "$O/truth_sample.txt"
  bcftools reheader -s "$O/truth_sample.txt" "$O/phased.vcf.gz" > "$O/phased.renamed.vcf.gz" && mv "$O/phased.renamed.vcf.gz" "$O/phased.vcf.gz"
  tabix -f -p vcf "$O/phased.vcf.gz"
  whatshap compare --names truth,calls --only-snvs --tsv-pairwise "$O/compare.tsv" "$TV" "$O/phased.vcf.gz" > "$O/compare.txt" 2>&1
  whatshap stats --chromosome chr20 --tsv "$O/stats.tsv" "$O/phased.vcf.gz" > "$O/stats.txt" 2>&1
  longphase haplotag -s "$O/phased.vcf.gz" -b "$A/aln.bam" -r "$REF" -t 16 -o "$O/haplotagged" > "$O/haplotag.log" 2>&1
  tagged=$(samtools view -c -F 0x904 -d HP:1 "$O/haplotagged.bam"); tagged=$((tagged + $(samtools view -c -F 0x904 -d HP:2 "$O/haplotagged.bam")))
  total=$(samtools view -c -F 0x904 "$O/haplotagged.bam")
  python3 - "$label" "$O/compare.tsv" "$O/stats.tsv" "$tagged" "$total" <<'PY' >> "$OUT/summary.tsv"
import csv, sys
label, cmp_, st, tagged, total = sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4]), int(sys.argv[5])
c = next(r for r in csv.DictReader(open(cmp_), delimiter="\t"))
s = next(r for r in csv.DictReader(open(st), delimiter="\t") if r["chromosome"] == "chr20")
n = int(c["all_assessed_pairs"]); sw = int(c["all_switches"])
print("\t".join(map(str, [label, c["covered_variants"], sw, c["all_switch_rate"], c["all_assessed_pairs"] and c.get("blockwise_hamming_rate", ""),
                          s["block_n50"], s["blocks"], f"{tagged / total:.4f}"])))
PY
  rm -f "$O/haplotagged.bam"*
done
column -t -s$'\t' "$OUT/summary.tsv"
