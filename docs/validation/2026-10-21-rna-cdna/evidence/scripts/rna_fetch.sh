#!/usr/bin/env bash
# RNA trial data: HG002 direct RNA (RNA004, SUP) and cDNA (PCS114, HAC) from ONT's UHRR_HG002_2026.06
# release (prefixes of the unaligned BAMs), and the GENCODE v49 annotation ONT's analysis used.
set -euo pipefail
R=${BENCH}; B=$R/bench; cd $R/data/rna
export PATH=$HOME/micromamba/envs/vareval/bin:$PATH
U=s3://ont-open-data/UHRR_HG002_2026.06/basecalls
get() { # name s3 bytes
  [[ -s $1.fq.gz ]] && return
  python3 $B/s3_parallel_cat.py $2 16 32 $3 > $1.prefix.bam.part && mv $1.prefix.bam.part $1.prefix.bam
  (samtools fastq -T MM,ML $1.prefix.bam 2>$1.fastq.err || true) | gzip -1 > $1.fq.gz.part
  mv $1.fq.gz.part $1.fq.gz; rm $1.prefix.bam
}
get drna_hg002 $U/dRNA/SUP/HG002/RNA004/HG002_RNA004_PolyA_1/HG002_RNA004_PolyA_1.bam 2e9
get cdna_hg002 $U/cDNA/HAC/HG002/PCS114/HG002_PCS114_PolyA_1/HG002_PCS114_PolyA_1.bam 2e9
[[ -s gencode.v49.annotation.gff3.gz ]] || curl -sSfL -o gencode.v49.annotation.gff3.gz \
  https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/gencode.v49.annotation.gff3.gz
echo "rna fetch done"
