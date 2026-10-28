#!/usr/bin/env bash
# chr20 HiFi test data: PacBio Revio 2022Q4 HG002/HG003/HG004 (reads whose primary alignment is on
# chr20 in PacBio's GRCh38 BAMs, streamed by region, methylation tags kept), GRCh38 no-alt, GIAB
# v4.2.1 small-variant truth (all three) and v5.0q SV truth (HG002), chr20 only.
set -euo pipefail
H=${BENCH}/hifi-data
export PATH=$HOME/micromamba/envs/vareval/bin:$PATH
P=https://downloads.pacbcloud.com/public/revio/2022Q4
G=https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio
cd $H/ref
if [[ ! -s GRCh38_no_alt.fa.fai ]]; then
  curl -sSfL https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz | gzip -dc > GRCh38_no_alt.fa
  samtools faidx GRCh38_no_alt.fa
fi
cd $H/reads
fetch() { # sample path
  [[ -s $1.chr20.fq.gz ]] && return
  samtools view -u -F 0x900 "$P/$2" chr20 | samtools fastq -T MM,ML - 2> $1.fastq.err | gzip -1 > $1.chr20.fq.gz.part
  mv $1.chr20.fq.gz.part $1.chr20.fq.gz
}
fetch HG002 HG002-rep1/analysis/HG002.m84011_220902_175841_s1.GRCh38.bam &
fetch HG003 HG003-rep1/analysis/HG003.m84010_220919_235306_s2.GRCh38.bam &
fetch HG004 HG004-rep1/analysis/HG004.m84010_220919_232145_s1.GRCh38.bam &
cd $H/truth
for s in HG002_NA24385_son HG003_NA24149_father HG004_NA24143_mother; do
  id=${s%%_*}; u=$G/$s/NISTv4.2.1/GRCh38/${id}_GRCh38_1_22_v4.2.1_benchmark
  [[ -s ${id}_chr20_v4.2.1.vcf.gz ]] || { bcftools view -r chr20 -Oz -o ${id}_chr20_v4.2.1.vcf.gz "$u.vcf.gz" && tabix -p vcf ${id}_chr20_v4.2.1.vcf.gz; }
  [[ -s ${id}_chr20_v4.2.1.bed ]] || curl -sSfL "${u}_noinconsistent.bed" | awk '$1=="chr20"' > ${id}_chr20_v4.2.1.bed
done
u=$G/HG002_NA24385_son/v5.0q/HG002_GRCh38_v5.0q_stvar
[[ -s HG002_chr20_v5.0q_stvar.vcf.gz ]] || { bcftools view -r chr20 -Oz -o HG002_chr20_v5.0q_stvar.vcf.gz "$u.vcf.gz" && tabix -p vcf HG002_chr20_v5.0q_stvar.vcf.gz; }
[[ -s HG002_chr20_v5.0q_stvar.bed ]] || curl -sSfL "$u.benchmark.bed" | awk '$1=="chr20"' > HG002_chr20_v5.0q_stvar.bed
wait
ls -la $H/reads $H/truth $H/ref
echo FETCH DONE
