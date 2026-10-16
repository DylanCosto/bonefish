#!/usr/bin/env bash
# HG002 genome-wide calls (bench/wgs_eval.sh: Clair3 + Sniffles2 on bonefish and minimap2 BAMs)
# scored against two GIAB benchmarks of hard regions:
#   CMRG v1.00 (273 challenging medically relevant genes): small variants with rtg vcfeval
#     (SNPs and indels separately), SVs with Truvari (--pctseq 0 --typeignore --dup-to-ins, as
#     the CMRG README advises: its duplications are encoded as insertions)
#   GIAB tandem-repeat benchmark v1.0.1: Clair3 PASS calls merged with Sniffles2 PASS SVs,
#     Truvari bench --sizemin 5 --pick ac on the Tier1 regions, then truvari refine (as the
#     benchmark README describes)
# Usage: bench/cmrg_tr_eval.sh OUTDIR     Env: TR_ONLY=1 (skip CMRG, keep its rows),
#        WGS_DIR (sample_eval.sh output, data/wgs/HG002), ALIGNERS ("minimap2 bonefish")
set -uo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
export PATH="${MAMBA_ROOT_PREFIX:-$HOME/micromamba}/envs/vareval/bin:$PATH"
OUT="${1:?outdir}"; mkdir -p "$OUT"
REF="$ROOT/data/ref/GRCh38_no_alt.fa"; SDF="$ROOT/data/ref/GRCh38_no_alt.sdf"
C="$ROOT/data/truth/cmrg"; TR="$ROOT/data/truth/tr"; W="${WGS_DIR:-$ROOT/data/wgs/HG002}"
CV="$C/HG002_GRCh38_CMRG_smallvar_v1.00.vcf.gz"; CB="$C/HG002_GRCh38_CMRG_smallvar_v1.00.bed"
for kind in snps indels; do
  [[ -s "$C/cmrg_$kind.vcf.gz.tbi" ]] || { bcftools view -v "$kind" -Oz -o "$C/cmrg_$kind.vcf.gz" "$CV"; tabix -p vcf "$C/cmrg_$kind.vcf.gz"; }
done
[[ -n "${TR_ONLY-}" ]] || echo -e "aligner\tbenchmark\ttype\ttp\tfp\tfn\tprecision\trecall\tf1" > "$OUT/summary.tsv"
for a in ${ALIGNERS:-minimap2 bonefish}; do
  O="$OUT/$a"; mkdir -p "$O"
  if [[ -z "${TR_ONLY-}" ]]; then
  # small variants in CMRG genes
  for kind in snps indels; do
    rm -rf "$O/cmrg_vcfeval_$kind"
    rtg vcfeval -b "$C/cmrg_$kind.vcf.gz" -c "$W/$a/calls_$kind.vcf.gz" -e "$CB" -t "$SDF" \
      -o "$O/cmrg_vcfeval_$kind" --threads 8 > /dev/null
    awk -v a="$a" -v k="$kind" '$1=="None"{printf "%s\tCMRG small\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", a, k, $2, $4, $5, $6, $7, $8}' \
      "$O/cmrg_vcfeval_$kind/summary.txt" >> "$OUT/summary.tsv"
  done
  # SVs in CMRG genes
  rm -rf "$O/cmrg_truvari"
  truvari bench -b "$C/HG002_GRCh38_CMRG_SV_v1.00.vcf.gz" -c "$W/$a/sv.vcf.gz" --includebed "$C/HG002_GRCh38_CMRG_SV_v1.00.bed" \
    --passonly --pctseq 0 --typeignore --dup-to-ins -o "$O/cmrg_truvari" > "$O/cmrg_truvari.log" 2>&1
  python3 -c "import json,sys; d=json.load(open(sys.argv[1])); print('\t'.join(map(str,[sys.argv[2],'CMRG SV','SV',d['TP-base'],d['FP'],d['FN'],round(d['precision'],4),round(d['recall'],4),round(d['f1'],4)])))" \
    "$O/cmrg_truvari/summary.json" "$a" >> "$OUT/summary.tsv"
  fi
  # tandem repeats: one call set of every size (Clair3 PASS + Sniffles2 PASS), sample renamed
  M="$O/merged_calls.vcf.gz"
  if [[ ! -s "$M.tbi" ]]; then
    echo SAMPLE > "$O/sample.txt"
    bcftools view -f PASS "$W/$a/clair3/merge_output.vcf.gz" | bcftools reheader -s "$O/sample.txt" | bcftools view -Oz -o "$O/clair3.pass.vcf.gz"
    bcftools view -f PASS "$W/$a/sv.vcf.gz" | bcftools reheader -s "$O/sample.txt" | bcftools view -Oz -o "$O/sv.pass.vcf.gz"
    tabix -f -p vcf "$O/clair3.pass.vcf.gz"; tabix -f -p vcf "$O/sv.pass.vcf.gz"
    # Truvari compares biallelic records only: split multi-allelic sites
    bcftools concat -a "$O/clair3.pass.vcf.gz" "$O/sv.pass.vcf.gz" | bcftools norm -m -any | bcftools sort -Oz -o "$M"
    tabix -p vcf "$M"; rm -f "$O"/clair3.pass.vcf.gz* "$O"/sv.pass.vcf.gz*
  fi
  rm -rf "$O/tr_bench"
  truvari bench -b "$TR/HG002_GRCh38_TandemRepeats_v1.0.1.vcf.gz" -c "$M" --includebed "$TR/HG002_GRCh38_TandemRepeats_v1.0.1_Tier1.bed.gz" \
    --sizemin 5 --pick ac -o "$O/tr_bench" > "$O/tr_bench.log" 2>&1
  python3 -c "import json,sys; d=json.load(open(sys.argv[1])); print('\t'.join(map(str,[sys.argv[2],'GIAB TR Tier1 (bench)','>=5 bp',d['TP-base'],d['FP'],d['FN'],round(d['precision'],4),round(d['recall'],4),round(d['f1'],4)])))" \
    "$O/tr_bench/summary.json" "$a" >> "$OUT/summary.tsv"
done
column -t -s$'\t' "$OUT/summary.tsv"
# refine (realigns candidate regions with MAFFT; slow genome-wide, so last)
for a in ${ALIGNERS:-minimap2 bonefish}; do
  O="$OUT/$a"
  truvari refine --use-original-vcfs --reference "$REF" --threads 8 "$O/tr_bench" > "$O/tr_refine.log" 2>&1
  python3 -c "import json,sys; d=json.load(open(sys.argv[1])); print('\t'.join(map(str,[sys.argv[2],'GIAB TR Tier1 (refine)','>=5 bp',d['TP-base'],d['FP'],d['FN'],round(d['precision'],4),round(d['recall'],4),round(d['f1'],4)])))" \
    "$O/tr_bench/refine.variant_summary.json" "$a" >> "$OUT/summary.tsv"
done
column -t -s$'\t' "$OUT/summary.tsv"
