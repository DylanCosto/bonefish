#!/usr/bin/env bash
# Score wf-human-variation outputs for the minimap2 and bonefish BAMs (bench/wfhv_compare.sh), HG002 chr20:
# small variants (rtg vcfeval vs GIAB v4.2.1), SVs (Truvari vs GIAB v5.0q), methylation (per-CpG agreement
# between the two runs), STRs (Straglr genotypes compared between runs), CNV calls (listed).
set -uo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
W="$ROOT/results/wfhv"; T="$ROOT/data/truth"; SDF="$ROOT/data/ref/GRCh38_no_alt.sdf"
export PATH="$HOME/micromamba/envs/vareval/bin:$PATH"
OUT="$W/score"; mkdir -p "$OUT"
TB="$T/HG002_chr20_v4.2.1.bed"; SVTV="$T/HG002_chr20_v5.0q_stvar.vcf.gz"; SVTB="$T/HG002_chr20_v5.0q_stvar.bed"
echo -e "aligner\ttype\ttruth_tp\tfp\tfn\tprecision\trecall\tf1" > "$OUT/small_variants.tsv"
echo -e "aligner\tTP\tFP\tFN\tprecision\trecall\tf1" > "$OUT/sv.tsv"
for a in minimap2 bonefish; do
  d="$W/$a"; s=HG002_$a
  for kind in snps indels; do
    bcftools view -f PASS -v $kind -Oz -o "$OUT/$a.$kind.vcf.gz" "$d/$s.wf_snp.vcf.gz" && tabix -f -p vcf "$OUT/$a.$kind.vcf.gz"
    rm -rf "$OUT/vcfeval_${a}_$kind"
    rtg vcfeval -b "$T/HG002_chr20_$kind.vcf.gz" -c "$OUT/$a.$kind.vcf.gz" -e "$TB" -t "$SDF" -o "$OUT/vcfeval_${a}_$kind" --threads 16 > /dev/null
    awk -v a=$a -v k=$kind '$1=="None"{print a"\t"k"\t"$2"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8}' "$OUT/vcfeval_${a}_$kind/summary.txt" >> "$OUT/small_variants.tsv"
  done
  rm -rf "$OUT/truvari_$a"
  truvari bench -b "$SVTV" -c "$d/$s.wf_sv.vcf.gz" --includebed "$SVTB" --passonly --pctseq 0 -o "$OUT/truvari_$a" > "$OUT/truvari_$a.log" 2>&1
  python3 - "$a" "$OUT/truvari_$a/summary.json" >> "$OUT/sv.tsv" <<'PY'
import json, sys
s = json.load(open(sys.argv[2]))
print(f"{sys.argv[1]}\t{s['TP-base']}\t{s['FP']}\t{s['FN']}\t{s['precision']:.4f}\t{s['recall']:.4f}\t{s['f1']:.4f}")
PY
done
# methylation: the workflow's bedMethyl (all CpGs, not phased)
mb() { ls "$W/$1"/*wf_mods*bedmethyl* 2>/dev/null | grep -v haplotype | head -1; }
A=$(mb minimap2); B=$(mb bonefish)
if [[ -n "$A" && -n "$B" ]]; then
  zcat -f "$A" > "$OUT/minimap2.bedmethyl"; zcat -f "$B" > "$OUT/bonefish.bedmethyl"
  micromamba run -n ontaln python "$ROOT/bench/modkit_compare.py" "$OUT/minimap2.bedmethyl" "$OUT/bonefish.bedmethyl" minimap2 bonefish 10 > "$OUT/mods.tsv"
fi
# STRs: Straglr genotypes per locus (allele copy numbers, read counts dropped), compared between the runs
python3 - "$W/minimap2/HG002_minimap2.wf_str.straglr.tsv" "$W/bonefish/HG002_bonefish.wf_str.straglr.tsv" > "$OUT/str.tsv" <<'PY'
import sys, csv, re
def load(p):
    d = {}
    for r in csv.reader(open(p), delimiter="\t"):
        if not r or r[0].startswith("#"): continue
        alleles = sorted(float(x) for x in re.findall(r"([0-9.]+)\(", r[4]))
        d[(r[0], r[1], r[2])] = alleles
    return d
a, b = load(sys.argv[1]), load(sys.argv[2])
common = [k for k in a if k in b]
same_n = [k for k in common if len(a[k]) == len(b[k])]
exact = sum(1 for k in same_n if a[k] == b[k])
within1 = sum(1 for k in same_n if all(abs(x - y) <= 1.0 for x, y in zip(a[k], b[k])))
print("loci_minimap2\tloci_bonefish\tcommon\tsame_allele_count\tidentical_alleles\talleles_within_1_unit")
print(f"{len(a)}\t{len(b)}\t{len(common)}\t{len(same_n)}\t{exact}\t{within1}")
for k in common:
    if len(a[k]) != len(b[k]) or any(abs(x - y) > 1.0 for x, y in zip(a[k], b[k])):
        print(f"# differs {k[0]}:{k[1]}-{k[2]} minimap2={a[k]} bonefish={b[k]}")
PY
for a in minimap2 bonefish; do
  echo "== $a"; zcat -f "$W/$a/HG002_$a.wf_cnv.vcf.gz" 2>/dev/null | grep -v '^#' | awk -F"\t" "{split(\$8,f,\";\"); s=\"\"; for(i in f) if (f[i] ~ /^(END|SVLEN|SVTYPE|CN)=/) s=s f[i] \";\"; print \$1, \$2, \$5, \$7, s}"
done > "$OUT/cnv.txt"
column -t "$OUT/small_variants.tsv"; column -t "$OUT/sv.tsv"; cat "$OUT/mods.tsv" 2>/dev/null | column -t; column -t "$OUT/str.tsv"; cat "$OUT/cnv.txt"
