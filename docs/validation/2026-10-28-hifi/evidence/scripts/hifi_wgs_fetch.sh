#!/usr/bin/env bash
# Whole-genome HiFi test data: HG002 Revio 2022Q4 rep1 HiFi reads (unaligned BAM, 45.8 GB) -> FASTQ.gz
# with methylation tags; GIAB v4.2.1 small-variant and v5.0q SV truth (whole genome).
set -euo pipefail
R=${BENCH}; H=$R/hifi-data; mkdir -p $H/wgs $H/truth/wgs
export PATH=$HOME/micromamba/envs/vareval/bin:$PATH
P=https://downloads.pacbcloud.com/public/revio/2022Q4/HG002-rep1/m84011_220902_175841_s1.hifi_reads.bam
G=https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son
cd $H/truth/wgs
for f in NISTv4.2.1/GRCh38/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz NISTv4.2.1/GRCh38/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi \
         NISTv4.2.1/GRCh38/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed \
         v5.0q/HG002_GRCh38_v5.0q_stvar.vcf.gz v5.0q/HG002_GRCh38_v5.0q_stvar.vcf.gz.tbi v5.0q/HG002_GRCh38_v5.0q_stvar.benchmark.bed; do
  [[ -s $(basename $f) ]] || curl -sSfL -o $(basename $f) $G/$f
done
cd $H/wgs
if [[ ! -s HG002.hifi_reads.bam ]]; then
  python3 $R/bench/s3_parallel_cat.py $P 16 32 > HG002.hifi_reads.bam.part
  mv HG002.hifi_reads.bam.part HG002.hifi_reads.bam
fi
[[ $(stat -c %s HG002.hifi_reads.bam) == 45756630027 ]] || { echo "size mismatch"; exit 1; }
if [[ ! -s HG002.hifi.fq.gz ]]; then
  samtools fastq -@ 4 -T MM,ML HG002.hifi_reads.bam 2> fastq.err | pigz -p 8 -1 > HG002.hifi.fq.gz.part
  mv HG002.hifi.fq.gz.part HG002.hifi.fq.gz
fi
tail -2 fastq.err
ls -la $H/wgs $H/truth/wgs
echo WGS FETCH DONE
