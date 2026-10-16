#!/usr/bin/env bash
# Whole-genome validation of bonefish against minimap2 on one GIAB ONT sample (R10.4.1, dorado
# sup v5.0.0, from s3://ont-open-data/giab_2025.01). Both aligners get the same reads (FASTQ
# streamed out of ONT's BAM: primary + unmapped records, no base-modification tags) and produce
# a coordinate-sorted, indexed BAM:
#   minimap2 -ax lr:hq | samtools sort --write-index      bonefish -a --sort -o aln.bam
# Then for each BAM: flagstat + bench/bam_summary.py (counts, identity, MAPQ, primary table),
# Clair3 (sup v5.0.0 model, longphase intermediate phasing), rtg vcfeval vs GIAB v4.2.1 overall
# and in GIAB v3.3 stratifications, Sniffles2 (+ Truvari vs GIAB v5.0q where a truth exists).
# Finally the primary alignments of the two aligners are compared and a report is written.
#
# Every step records wall / CPU / peak RSS (GNU time), and bench/monitor.sh samples CPU, memory,
# swap, iowait and disk every 10 s into OUT/resources.csv with the current stage.
# Stages are resumable (finished outputs are skipped). BAMs and reads are deleted when no longer
# needed (KEEP_BAM=1 / KEEP_READS=1 keep them).
#
# Usage: bench/sample_eval.sh SAMPLE      (HG005 | HG002)
# Env: THREADS (16), SMOKE=N (pipeline test on the first N reads, in data/wgs/SAMPLE_smoke),
#      SUFFIX (output in data/wgs/SAMPLE$SUFFIX), NEED_ALIGN_GB (free space required to align, 150)
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
SAMPLE="${1:?sample}"
THREADS="${THREADS:-16}"
VE="${MAMBA_ROOT_PREFIX:-$HOME/micromamba}/envs/vareval/bin"
export PATH="$VE:${MAMBA_ROOT_PREFIX:-$HOME/micromamba}/envs/ontaln/bin:$PATH"
MM2="${MAMBA_ROOT_PREFIX:-$HOME/micromamba}/envs/ontaln/bin/minimap2"
BONEFISH="${BONEFISH:-$ROOT/bonefish/target/release/bonefish}"
REF="$ROOT/data/ref/GRCh38_no_alt.fa"; SDF="$ROOT/data/ref/GRCh38_no_alt.sdf"
MMI="$ROOT/data/index/mm2.lrhq.mmi"; RAW="$ROOT/data/index/bonefish.lrhq.raw"
MODEL="$VE/models/r1041_e82_400bps_sup_v500"
T="$ROOT/data/truth"; ST="$T/strat"
S3=s3://ont-open-data/giab_2025.01/basecalling/sup
STRATA=(alldifficultregions notinalldifficultregions segdups lowmappabilityall AllTandemRepeatsandHomopolymers_slop5)

FRAC=""; SVT=""; SVB=""
case "$SAMPLE" in
  HG005) SRC=$S3/HG005/PAW87816/calls.sorted.bam
         GIAB=https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/ChineseTrio/HG005_NA24631_son/NISTv4.2.1/GRCh38
         TV="$T/HG005_GRCh38_1_22_v4.2.1_benchmark.vcf.gz"; TB="$T/HG005_GRCh38_1_22_v4.2.1_benchmark.bed" ;;
  HG002) SRC=$S3/HG002/PAW71238/calls.sorted.bam
         FRAC=0.6   # ~50x flowcell -> ~30x
         TV="$T/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz"; TB="$T/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed"
         SVT="$T/HG002_GRCh38_v5.0q_stvar.vcf.gz"; SVB="$T/HG002_GRCh38_v5.0q_stvar.benchmark.bed" ;;
  *) echo "unknown sample $SAMPLE"; exit 1 ;;
esac

O="$ROOT/data/wgs/$SAMPLE${SUFFIX-}${SMOKE:+_smoke}"; mkdir -p "$O"
O="$(cd "$O" && pwd)"   # absolute (Clair3 breaks on the volume path's space otherwise)
READS="$O/reads.fq.gz"
log() { echo "[$(date '+%F %T')] $SAMPLE: $*" | tee -a "$O/run.log"; }
stage() { echo "$1" > "$O/stage"; log "stage: $1"; }
TFMT="%e %U %S %M"   # wall s, user s, sys s, peak RSS KB
need_disk() { local free; free=$(df -BG --output=avail "$O/" | tail -1 | tr -dc '0-9')
  (( free >= $1 )) || { log "ABORT: only ${free} GB free, need $1 GB for $2"; exit 1; }; }

MON_CSV="$O/resources.csv" MON_STAGE="$O/stage" "$ROOT/bench/monitor.sh" 10 1000000 >> "$O/monitor_alerts.log" 2>&1 &
MON_PID=$!
trap 'kill $MON_PID 2>/dev/null || true' EXIT
log "start (bonefish $("$BONEFISH" --version | cut -d" " -f2), $("$MM2" --version | sed 's/^/minimap2 /'))"

# --- truth ------------------------------------------------------------------------------------
stage truth
if [[ ! -s "$TV.tbi" ]]; then
  for f in "$(basename "$TV")" "$(basename "$TV").tbi" "$(basename "$TB")"; do curl -fsSL --retry 3 -o "$T/$f" "$GIAB/$f"; done
fi
for kind in snps indels; do
  [[ -s "$O/truth_$kind.vcf.gz.tbi" ]] || { bcftools view -v "$kind" -Oz -o "$O/truth_$kind.vcf.gz" "$TV"; tabix -p vcf "$O/truth_$kind.vcf.gz"; }
done
for s in "${STRATA[@]}"; do
  [[ -s "$O/bench_$s.bed" ]] || "$ROOT/bench/bed_intersect.py" "$TB" "$ST/GRCh38_$s.bed.gz" > "$O/bench_$s.bed"
done

# --- reads: stream ONT's BAM from S3 into a bgzipped FASTQ --------------------------------------
if [[ ! -s "$READS" ]]; then
  stage download
  need_disk 130 "reads"
  aws configure set default.s3.max_concurrent_requests 32 2>/dev/null || true
  sub=(); [[ -n "$FRAC" ]] && sub=(-s "7${FRAC#0}")
  if [[ -n "${SMOKE-}" ]]; then
    { aws s3 cp --no-sign-request --only-show-errors "$SRC" - 2>/dev/null | samtools view -u -F 0x900 - 2>/dev/null \
        | samtools fastq - 2>/dev/null | head -n $((4 * SMOKE)) | bgzip > "$READS.part"; } || true
    echo "0 0 0 0" > "$O/download.time"
  else
    # A parallel (multipart) download to a file is ~2x faster than streaming S3 to a pipe.
    need_disk 230 "download + FASTQ"
    [[ -s "$O/input.bam" ]] || { /usr/bin/time -f "$TFMT" -o "$O/download.time" \
      aws s3 cp --no-sign-request --only-show-errors "$SRC" "$O/input.bam.part" && mv "$O/input.bam.part" "$O/input.bam"; }
    log "downloaded $(du -h "$O/input.bam" | cut -f1) in $(cut -d' ' -f1 "$O/download.time") s"
    stage convert
    /usr/bin/time -f "$TFMT" -o "$O/convert.time" bash -c "set -o pipefail
      samtools view -u -@ 4 -F 0x900 ${sub[*]} '$O/input.bam' \
        | samtools fastq -@ 4 - 2>/dev/null \
        | bgzip -@ 8 > '$READS.part'"
    rm -f "$O/input.bam"
  fi
  mv "$READS.part" "$READS"
  log "reads: $(du -h "$READS" | cut -f1) in $(cut -d' ' -f1 "$O/download.time") s"
fi

# --- per aligner --------------------------------------------------------------------------------
for A in minimap2 bonefish; do
  D="$O/$A"; mkdir -p "$D"; BAM="$D/aln.bam"; CALLS="$D/clair3/merge_output.vcf.gz"
  if [[ ! -s "$D/bam_summary.tsv" && ! -s "$BAM.bai" ]]; then
    stage "align_$A"
    need_disk "${NEED_ALIGN_GB:-150}" "align $A"
    rm -f "$BAM" "$BAM.bai" "$D"/sort.tmp* "$BAM".tmp*
    t0=$(date +%s.%N)
    if [[ $A == minimap2 ]]; then
      cat "$MMI" > /dev/null
      bash -c "set -o pipefail
        /usr/bin/time -f '$TFMT' -o '$D/align.time' '$MM2' -ax lr:hq -t $THREADS -K 50M '$MMI' '$READS' 2> '$D/align.log' \
          | /usr/bin/time -f '$TFMT' -o '$D/sort.time' samtools sort -@ $THREADS -m 384M --write-index -T '$D/sort.tmp' -o '$BAM##idx##$BAM.bai' -"
    else
      cat "$RAW" > /dev/null
      /usr/bin/time -f "$TFMT" -o "$D/align.time" "$BONEFISH" -a -t "$THREADS" -K 50M --sort -o "$BAM" "$RAW" "$READS" 2> "$D/align.log"
    fi
    awk -v a="$(date +%s.%N)" -v b="$t0" 'BEGIN{printf "%.1f\n", a - b}' > "$D/to_sorted_bam.seconds"
    log "$A: reads -> sorted BAM in $(cut -d. -f1 "$D/to_sorted_bam.seconds") s (aligner: $(cat "$D/align.time"))"
    ls -l "$BAM" | awk '{print $5}' > "$D/bam_bytes"
  fi
  if [[ ! -s "$D/bam_summary.tsv" ]]; then
    stage "bamstats_$A"
    samtools flagstat -@ 8 "$BAM" > "$D/flagstat.txt"
    /usr/bin/time -f "$TFMT" -o "$D/bamstats.time" bash -c "set -o pipefail
      samtools view -@ 4 '$BAM' | cut -f1-6,12- | python3 '$ROOT/bench/bam_summary.py' '$D'"
    log "$A: $(tr '\n' ' ' < "$D/bam_summary.tsv" | cut -c1-400)"
  fi
  if [[ ! -s "$CALLS" ]]; then
    stage "clair3_$A"
    rm -rf "$D/clair3"
    /usr/bin/time -f "$TFMT" -o "$D/clair3.time" run_clair3.sh --bam_fn="$BAM" --ref_fn="$REF" \
      --threads="$THREADS" --platform=ont --model_path="$MODEL" --output="$D/clair3" \
      --use_longphase_for_intermediate_phasing=True --longphase="$VE/longphase" > "$D/clair3.log" 2>&1
    log "$A: clair3 $(cat "$D/clair3.time")"
  fi
  if [[ ! -s "$D/strat.tsv" ]]; then
    stage "vcfeval_$A"
    {
      echo -e "set\tvalue"
      echo -e "calls_total\t$(bcftools view -H "$CALLS" | wc -l)"
      echo -e "calls_pass_snps\t$(bcftools view -H -f PASS -v snps "$CALLS" | wc -l)"
      echo -e "calls_pass_indels\t$(bcftools view -H -f PASS -v indels "$CALLS" | wc -l)"
    } > "$D/clair3_counts.tsv"
    echo -e "stratum\ttype\ttp\tfp\tfn\tprecision\trecall\tf1" > "$D/accuracy.tsv"
    for kind in snps indels; do
      bcftools view -v "$kind" -f PASS -Oz -o "$D/calls_$kind.vcf.gz" "$CALLS"; tabix -f -p vcf "$D/calls_$kind.vcf.gz"
      rm -rf "$D/vcfeval_$kind"
      rtg vcfeval -b "$O/truth_$kind.vcf.gz" -c "$D/calls_$kind.vcf.gz" -e "$TB" -t "$SDF" -o "$D/vcfeval_$kind" --threads "$THREADS" > /dev/null
      awk -v k="$kind" '$1=="None"{printf "all\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", k, $2, $4, $5, $6, $7, $8}' "$D/vcfeval_$kind/summary.txt" >> "$D/accuracy.tsv"
    done
    log "$A accuracy: $(tail -n +2 "$D/accuracy.tsv" | awk '{printf "%s F1 %s (FP %s FN %s); ", $2, $8, $4, $5}')"
    echo -e "stratum\ttype\ttp\tfp\tfn\tprecision\trecall\tf1" > "$D/strat.tsv.part"
    for s in "${STRATA[@]}"; do
      for kind in snps indels; do
        d="$D/strat_${s}_$kind"; rm -rf "$d"
        rtg vcfeval -b "$O/truth_$kind.vcf.gz" -c "$D/calls_$kind.vcf.gz" -e "$O/bench_$s.bed" -t "$SDF" -o "$d" --threads "$THREADS" > /dev/null
        awk -v s="$s" -v k="$kind" '$1=="None"{printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", s, k, $2, $4, $5, $6, $7, $8}' "$d/summary.txt" >> "$D/strat.tsv.part"
        rm -rf "$d"
      done
    done
    mv "$D/strat.tsv.part" "$D/strat.tsv"
  fi
  if [[ ! -s "$D/sv.vcf.gz.tbi" ]]; then
    stage "sniffles_$A"
    /usr/bin/time -f "$TFMT" -o "$D/sniffles.time" sniffles --input "$BAM" --vcf "$D/sv.vcf" --reference "$REF" \
      --threads "$THREADS" > "$D/sniffles.log" 2>&1
    bcftools sort -Oz -o "$D/sv.vcf.gz" "$D/sv.vcf" 2>/dev/null && tabix -f -p vcf "$D/sv.vcf.gz" && rm -f "$D/sv.vcf"
    log "$A: sniffles $(cat "$D/sniffles.time"), $(bcftools view -H -f PASS "$D/sv.vcf.gz" | wc -l) PASS SVs"
  fi
  if [[ -n "$SVT" && ! -s "$D/sv_accuracy.tsv" ]]; then
    stage "truvari_$A"
    echo -e "matching\ttp\tfp\tfn\tprecision\trecall\tf1" > "$D/sv_accuracy.tsv"
    for mode in strict multi; do
      extra=(); [[ $mode == multi ]] && extra=(--pick multi)
      rm -rf "$D/truvari_$mode"
      truvari bench -b "$SVT" -c "$D/sv.vcf.gz" --includebed "$SVB" --passonly --pctseq 0 "${extra[@]}" \
        -o "$D/truvari_$mode" > "$D/truvari_$mode.log" 2>&1
      python3 -c "import json,sys; s=json.load(open('$D/truvari_$mode/summary.json')); print('\t'.join(['$mode',str(s['TP-base']),str(s['FP']),str(s['FN'])]+[f\"{s[k]:.4f}\" for k in ('precision','recall','f1')]))" >> "$D/sv_accuracy.tsv"
    done
    log "$A SV: $(tail -n +2 "$D/sv_accuracy.tsv" | awk '{printf "%s F1 %s (TP %s FP %s); ", $1, $7, $2, $3}')"
  fi
  if [[ -z "${KEEP_BAM-}" && -s "$D/strat.tsv" && -s "$D/sv.vcf.gz.tbi" ]]; then
    rm -f "$BAM" "$BAM.bai"; rm -rf "$D/clair3/tmp"
  fi
done

# --- aligner comparison and report ------------------------------------------------------------
stage compare
[[ -s "$O/compare_primary.tsv" ]] || python3 "$ROOT/bench/compare_tables.py" "$O/minimap2/primary.tsv.gz" "$O/bonefish/primary.tsv.gz" > "$O/compare_primary.tsv"
if [[ ! -s "$O/sv_agreement.tsv" ]]; then
  rm -rf "$O/truvari_agreement"
  truvari bench -b "$O/minimap2/sv.vcf.gz" -c "$O/bonefish/sv.vcf.gz" --passonly --pctseq 0 -o "$O/truvari_agreement" > /dev/null 2>&1
  python3 -c "import json; s=json.load(open('$O/truvari_agreement/summary.json')); print('minimap2_svs_found_by_bonefish\t%d\nbonefish_only\t%d\nminimap2_only\t%d\nagreement_f1\t%.4f' % (s['TP-base'], s['FP'], s['FN'], s['f1']))" > "$O/sv_agreement.tsv"
fi
stage report
python3 "$ROOT/bench/sample_report.py" "$O" > "$O/report.md"
log "report: $O/report.md"
[[ -n "${KEEP_READS-}" ]] || rm -f "$READS"
stage done
log "done"
