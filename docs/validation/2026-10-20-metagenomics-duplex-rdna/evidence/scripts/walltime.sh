#!/usr/bin/env bash
# Wall time without the single-threaded `samtools view -b` on the output (SAM to /dev/null).
R=${BENCH}; M=$R/data/meta; export PATH=$HOME/micromamba/envs/ontaln/bin:$PATH
D=$R/data/duplex/11_15_22_R1041_Duplex_HG002_3_Dorado_v0.1.1_400bps_sup_stereo_duplex_pass.fastq.gz
t() { /usr/bin/time -f "$1: %e s wall, %U+%S CPU-s" -a -o walltime.txt "${@:2}" > /dev/null 2>/dev/null; }
: > walltime.txt
t "meta bonefish" $R/bin/bonefish-v017 -a -y -t 16 $M/meta.bonefish.raw $M/fecal_reads.fq.gz
t "meta minimap2" minimap2 -ax lr:hq -y -t 16 $M/meta.mm2.mmi $M/fecal_reads.fq.gz
t "duplex bonefish" $R/bin/bonefish-v017 -a -y -t 16 $R/data/index/bonefish.lrhq.raw $D
t "gunzip -c fecal reads" gzip -dc $M/fecal_reads.fq.gz
t "gunzip -c duplex reads" gzip -dc $D
echo DONE >> walltime.txt
