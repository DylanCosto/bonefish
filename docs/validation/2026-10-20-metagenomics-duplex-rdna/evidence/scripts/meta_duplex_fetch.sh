#!/usr/bin/env bash
# Downloads for the metagenomics (Zymo fecal reference, ONT open data) and duplex (HG002, HPRC/UCSC) checks.
set -euo pipefail
R=${BENCH}; B=$R/bench
export PATH=$HOME/micromamba/envs/vareval/bin:$PATH
cd $R/data/meta
Z=zymo_fecal_2025.05
[[ -s contigs.fa.gz ]] || curl -sSf -o contigs.fa.gz https://ont-open-data.s3.amazonaws.com/$Z/analysis/EigthFC/bgzip_contigs.fa.gz
for f in bgzip_contigs.fa.gz.circular.list dastool/das_DASTool_contig2bin.tsv dastool/das_DASTool_summary.tsv; do
  [[ -s $(basename $f) ]] || curl -sSf -o $(basename $f) https://ont-open-data.s3.amazonaws.com/$Z/analysis/EigthFC/$f; done
curl -sSf -o Analysis_commands.sh https://ont-open-data.s3.amazonaws.com/$Z/analysis/Analysis_commands.sh
if [[ ! -s fecal_reads.fq.gz ]]; then
  # first 12 GB of one flowcell's sup BAM (unaligned, with modified-base tags); the cut record at the end is dropped
  python3 $B/s3_parallel_cat.py s3://ont-open-data/$Z/basecalls/PAU85136/basecalls.sup_5.0.0.4mC5mC6mA.bam 16 32 12e9 > prefix.bam.part
  mv prefix.bam.part prefix.bam
  (samtools fastq -T MM,ML prefix.bam 2>fastq.err || true) | gzip -1 > fecal_reads.fq.gz.part
  mv fecal_reads.fq.gz.part fecal_reads.fq.gz; rm prefix.bam
fi
echo "meta done"
cd $R/data/duplex
D=submissions/0CB931D5-AE0C-4187-8BD8-B3A9C9BFDADE--UCSC_HG002_R1041_Duplex_Dorado/Dorado_v0.1.1/stereo_duplex
F=11_15_22_R1041_Duplex_HG002_3_Dorado_v0.1.1_400bps_sup_stereo_duplex_pass.fastq.gz
[[ -s $F ]] || { python3 $B/s3_parallel_cat.py s3://human-pangenomics/$D/$F 16 32 > $F.part && mv $F.part $F; }
echo "duplex done"
