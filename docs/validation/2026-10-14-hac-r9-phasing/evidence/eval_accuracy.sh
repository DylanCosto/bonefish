#!/usr/bin/env bash
# Phase 3 accuracy harness: align one chromosome's reads (default chr20, ~55x from flowcell
# PAW70337; CHROM=chr22 READS=... for others), call small variants with Clair3 (ONT R10.4.1 sup
# v5.0.0 model), and score against GIAB HG002 v4.2.1 on that chromosome with rtg vcfeval (SNPs
# and indels separately). Also compares primary alignments with a reference BAM when given,
# and runs the SV check (bench/eval_sv.sh).
#
# Usage: bench/eval_accuracy.sh LABEL "ALIGNER_ARGS" [REFERENCE_BAM]
#   e.g. bench/eval_accuracy.sh baseline ""
#        bench/eval_accuracy.sh zdrop200 "-z 200" results/accuracy/baseline/aln.bam
# Env: BONEFISH (binary, default bonefish/target/release/bonefish), PRESET (default dna),
#      ALIGNER=minimap2 (align with minimap2 -x lr:hq instead, the baseline), THREADS (default 16),
#      CHROM (default chr20), READS (default: the chr20 reads), KEEP_BAM, CLAIR3_MODEL (model dir name).
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
LABEL="${1:?label}"; ARGS="${2-}"; REFBAM="${3-}"
THREADS="${THREADS:-16}"
BIN="${BONEFISH:-$ROOT/bonefish/target/release/bonefish}"
VE="${MAMBA_ROOT_PREFIX:-$HOME/micromamba}/envs/vareval/bin"
export PATH="$VE:$PATH"

REF="$ROOT/data/ref/GRCh38_no_alt.fa"
IDX="${IDX:-$ROOT/data/index/bonefish.lrhq.raw}"
PRESET="${PRESET:-dna}"
if [[ "${ALIGNER-}" == minimap2 ]]; then
  BIN="${MAMBA_ROOT_PREFIX:-$HOME/micromamba}/envs/ontaln/bin/minimap2"
  IDX="$ROOT/data/index/mm2.lrhq.mmi"; PRESET=lr:hq
fi
CHROM="${CHROM:-chr20}"
READS="${READS:-$ROOT/data/reads/hg002_PAW70337_chr20.fq.gz}"
MODEL="$VE/models/${CLAIR3_MODEL:-r1041_e82_400bps_sup_v500}"   # PyTorch model bundled with Clair3 2.x
TRUTH_DIR="$ROOT/data/truth"
SDF="$ROOT/data/ref/GRCh38_no_alt.sdf"
OUT="$ROOT/results/accuracy/$LABEL"
mkdir -p "$OUT"
log() { echo "[$(date +%T)] $*" | tee -a "$OUT/eval.log"; }

# --- one-time inputs: per-chromosome truth subset, rtg reference SDF -----------------------
TV="$TRUTH_DIR/HG002_${CHROM}_v4.2.1.vcf.gz"; TB="$TRUTH_DIR/HG002_${CHROM}_v4.2.1.bed"
if [[ ! -s "$TV.tbi" ]]; then
  bcftools view -r "$CHROM" -Oz -o "$TV" "$TRUTH_DIR/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz"
  tabix -p vcf "$TV"
  awk -v c="$CHROM" '$1==c' "$TRUTH_DIR/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed" > "$TB"
fi
for kind in snps indels; do
  [[ -s "$TRUTH_DIR/HG002_${CHROM}_$kind.vcf.gz.tbi" ]] || {
    bcftools view -v "$kind" -Oz -o "$TRUTH_DIR/HG002_${CHROM}_$kind.vcf.gz" "$TV"
    tabix -p vcf "$TRUTH_DIR/HG002_${CHROM}_$kind.vcf.gz"; }
done
[[ -d "$SDF" ]] || rtg format -o "$SDF" "$REF" > /dev/null

# --- 1. align + sort (accuracy only; speed is measured separately with checksummed SAM) ---
BAM="$OUT/aln.bam"
if [[ ! -s "$BAM.bai" ]]; then
  log "align: $BIN -ax $PRESET $ARGS"
  cat "$IDX" > /dev/null
  set -o pipefail
  /usr/bin/time -f "%e s wall, %U+%S cpu, %M KB" -o "$OUT/align.time" \
    "$BIN" -ax "$PRESET" -t "$THREADS" -K 50M $ARGS "$IDX" "$READS" 2> "$OUT/align.log" \
    | samtools sort -@ 4 -m 1500M -o "$BAM" -
  samtools index -@ 4 "$BAM"
  log "align done: $(cat "$OUT/align.time")"
fi

# --- 2. Clair3 on the chromosome -------------------------------------------------------------------
CALLS="$OUT/clair3/merge_output.vcf.gz"
if [[ ! -s "$CALLS" ]]; then
  log "clair3"
  /usr/bin/time -f "%e s wall" -o "$OUT/clair3.time" \
    run_clair3.sh --bam_fn="$BAM" --ref_fn="$REF" --threads="$THREADS" --platform=ont \
      --model_path="$MODEL" --output="$OUT/clair3" --ctg_name="$CHROM" > "$OUT/clair3.log" 2>&1
  log "clair3 done: $(cat "$OUT/clair3.time")"
fi

# --- 3. rtg vcfeval against GIAB, SNPs and indels separately -------------------------------
summary="$OUT/accuracy.tsv"
echo -e "label\ttype\ttruth_tp\tfp\tfn\tprecision\trecall\tf1" > "$summary"
for kind in snps indels; do
  bcftools view -v "$kind" -f PASS -Oz -o "$OUT/calls_$kind.vcf.gz" "$CALLS"
  tabix -f -p vcf "$OUT/calls_$kind.vcf.gz"
  rm -rf "$OUT/vcfeval_$kind"
  rtg vcfeval -b "$TRUTH_DIR/HG002_${CHROM}_$kind.vcf.gz" -c "$OUT/calls_$kind.vcf.gz" -e "$TB" \
    -t "$SDF" -o "$OUT/vcfeval_$kind" --threads "$THREADS" > /dev/null
  # "None" row = all calls (no score threshold)
  awk -v l="$LABEL" -v k="$kind" '$1=="None"{printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", l, k, $2, $4, $5, $6, $7, $8}' \
    "$OUT/vcfeval_$kind/summary.txt" >> "$summary"
done
log "accuracy:"; column -t "$summary" | tee -a "$OUT/eval.log"

# --- 4. primary-alignment concordance with a reference BAM ---------------------------------
if [[ -n "$REFBAM" ]]; then
  python3 "$ROOT/bench/compare_primary.py" "$REFBAM" "$BAM" | tee "$OUT/concordance.txt" | tee -a "$OUT/eval.log"
fi

# --- 5. structural variants (Sniffles2 + Truvari vs GIAB v5.0q) ----------------------------
CHROM="$CHROM" "$ROOT/bench/eval_sv.sh" "$LABEL" | tee -a "$OUT/eval.log"

# Keep disk in check: drop this run's BAM once everything above has used it (KEEP_BAM=1 keeps it).
[[ -n "${KEEP_BAM-}" ]] || rm -f "$BAM" "$BAM.bai"
