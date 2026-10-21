#!/usr/bin/env bash
# DNA guard: rna-opt build vs v0.1.7 release on the genome-wide HG001 sample (output MD5, 3 alternating timing rounds)
D=${DATA}; IDX=$D/index/bonefish.lrhq.raw; R=$D/reads/hg001_wg_1x.fq.gz
B1=${BENCH}/bin/bonefish-v017; B2=${BENCH}/bin/bf-rna-opt1
cd $(dirname $0); : > report.txt
python3 -c "
import os
for f in ['$D/../data/rna/bf.rna.raw','$D/../data/rna/bf.cdna.raw']:
    if os.path.exists(f):
        fd=os.open(f,os.O_RDONLY); os.posix_fadvise(fd,0,0,os.POSIX_FADV_DONTNEED); os.close(fd)"
cat $IDX > /dev/null
for i in 1 2 3; do for b in $B1 $B2; do
  /usr/bin/time -f "%e %U %S" -o t.tmp $b -a -t 16 $IDX $R 2>/dev/null | grep -v '^@PG' | md5sum | cut -d' ' -f1 > m.tmp
  read w u s < t.tmp; echo "$i $(basename $b) wall=$w cpu=$(echo "$u+$s"|bc) md5=$(cat m.tmp)" >> report.txt
done; done
echo DONE >> report.txt
