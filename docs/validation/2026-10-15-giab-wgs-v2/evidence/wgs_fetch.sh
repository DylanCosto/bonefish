#!/usr/bin/env bash
# Stream one GIAB ONT sample (R10.4.1 sup, s3://ont-open-data/giab_2025.01) to FASTQ without storing
# the BAM: primary + unmapped records, optional subsample, as in sample_eval.sh's convert stage.
# Usage: [PARALLEL=N] bench/wgs_fetch.sh SAMPLE OUTDIR   (HG002: 60% of PAW71238, seed 7; HG005: all of PAW87816)
set -euo pipefail
SAMPLE=$1; O=$2
E=$HOME/micromamba/envs; export PATH=$E/vareval/bin:$E/ontaln/bin:$PATH
S3=s3://ont-open-data/giab_2025.01/basecalling/sup
case $SAMPLE in
  HG002) SRC=$S3/HG002/PAW71238/calls.sorted.bam; SUB=(-s 7.6) ;;
  HG005) SRC=$S3/HG005/PAW87816/calls.sorted.bam; SUB=() ;;
esac
# PARALLEL=N: N parallel ranged GETs (bench/s3_parallel_cat.py); a single stream tops out at ~10-15 MB/s
if [[ -n "${PARALLEL-}" ]]; then GET="python3 -I $(cd "$(dirname "$0")" && pwd)/s3_parallel_cat.py '$SRC' $PARALLEL 8"
else GET="aws s3 cp --no-sign-request --only-show-errors '$SRC' -"; fi
mkdir -p "$O"; [[ -s $O/reads.fq.gz ]] && { echo "already there"; exit 0; }
echo "[$(date '+%F %T')] start $SRC" >> "$O/fetch.log"
/usr/bin/time -f "%e %U %S %M" -o "$O/fetch.time" bash -c "set -o pipefail
  $GET \
    | samtools view -u -@ 4 -F 0x900 ${SUB[*]} - \
    | samtools fastq -@ 4 - 2>/dev/null | bgzip -@ 8 > '$O/reads.fq.gz.part'"
mv "$O/reads.fq.gz.part" "$O/reads.fq.gz"
echo "[$(date '+%F %T')] done $(du -h "$O/reads.fq.gz" | cut -f1), $(cut -d' ' -f1 "$O/fetch.time") s" >> "$O/fetch.log"
