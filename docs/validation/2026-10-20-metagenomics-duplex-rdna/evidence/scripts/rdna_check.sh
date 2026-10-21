#!/usr/bin/env bash
# rDNA reads: how bonefish and minimap2 place reads from the 45S rDNA arrays.
# GRCh38 holds one rDNA unit region on chrUn_GL000220v1; T2T-CHM13 (hs1) models the arrays on the
# acrocentric short arms. Reads with any alignment to GL000220 (either aligner) are the rDNA set.
set -euo pipefail
R=${BENCH}; O=$R/results/rdna; cd $O
export PATH=$HOME/micromamba/envs/vareval/bin:$HOME/micromamba/envs/ontaln/bin:$PATH
BF=${BF:?bonefish binary}; T=16
READS=$R/data/reads/hg001_wg_1x.fq.gz
[[ -s bf_grch38.bam ]] || $BF -a -t $T $R/data/index/bonefish.lrhq.raw $READS | samtools view -b -o bf_grch38.bam -
[[ -s mm_grch38.bam ]] || minimap2 -ax lr:hq -t $T $R/data/index/mm2.lrhq.mmi $READS 2>mm_grch38.err | samtools view -b -o mm_grch38.bam -
for a in bf mm; do samtools view ${a}_grch38.bam | awk '$3=="chrUn_GL000220v1"{print $1}' ; done | sort -u > rdna_names.txt
echo "rDNA read set: $(wc -l < rdna_names.txt) reads"
samtools view -N rdna_names.txt -b -o rdna_reads_from_bf.bam bf_grch38.bam
samtools fastq -F 0x900 rdna_reads_from_bf.bam 2>/dev/null | gzip > rdna_reads.fq.gz
[[ -s bf_hs1.bam ]] || $BF -a -t $T $R/data/index/extra/hs1.bonefish.raw rdna_reads.fq.gz | samtools view -b -o bf_hs1.bam -
[[ -s mm_hs1.bam ]] || minimap2 -ax lr:hq -t $T $R/data/ref/extra/hs1.fa rdna_reads.fq.gz 2>mm_hs1.err | samtools view -b -o mm_hs1.bam -
for ref in grch38 hs1; do for a in bf mm; do
  f=${a}_${ref}.bam; [[ $ref == grch38 ]] && f="-N rdna_names.txt $f"
  samtools view -F 0x900 $f | awk -v tag="$a $ref" '
    {n++; if($2%8>=4){u++; next}; c[$3]++; if($5>=20)q20++; if($5==0)q0++}
    END{printf "%s primary=%d unmapped=%d mapq0=%d mapq>=20=%d contigs:", tag, n, u, q0, q20;
        for(k in c) if(c[k]>=5) printf " %s=%d", k, c[k]; print ""}'
done; done | tee summary.txt
