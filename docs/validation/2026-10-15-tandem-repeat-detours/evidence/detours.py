#!/usr/bin/env python3
"""Count deletion-insertion detours in tandem repeats: primary alignments with an insertion and a
deletion of >= MIN bp each within WIN reference bp of each other (the pattern a misplaced anchor in a
repeat produces), inside the given BED regions (e.g. GIAB TR Tier 1).
Usage: detours.py BAM REGIONS.bed CHROM [MIN=30] [WIN=300]
"""
import bisect, re, subprocess, sys
bam, bed, chrom = sys.argv[1:4]
MIN = int(sys.argv[4]) if len(sys.argv) > 4 else 30
WIN = int(sys.argv[5]) if len(sys.argv) > 5 else 300
iv = sorted((int(l.split()[1]), int(l.split()[2])) for l in open(bed) if l.split()[0] == chrom)
starts = [a for a, _ in iv]
def in_tr(p):
    i = bisect.bisect_right(starts, p) - 1
    return i >= 0 and iv[i][0] <= p < iv[i][1]
p = subprocess.Popen(["samtools", "view", "-F", "0x904", bam, chrom], stdout=subprocess.PIPE, text=True)
reads = detour_reads = 0; loci = set()
for l in p.stdout:
    f = l.split("\t", 6)
    r = int(f[3]); ev = []
    for n, op in re.findall(r"(\d+)([MIDNSHP=X])", f[5]):
        n = int(n)
        if op in "M=XN": r += n
        elif op == "D":
            if n >= MIN: ev.append(("D", r, n))
            r += n
        elif op == "I" and n >= MIN: ev.append(("I", r, n))
    reads += 1
    hit = False
    for i in range(len(ev)):
        for j in range(i + 1, len(ev)):
            if ev[j][1] - ev[i][1] > WIN: break
            if ev[i][0] != ev[j][0] and in_tr(ev[i][1]) and in_tr(ev[j][1]):
                hit = True; loci.add(ev[i][1] // 500)
    detour_reads += hit
print(f"{bam.split('/')[-2]}\t{chrom}\treads {reads}\tdetour reads {detour_reads}\tloci {len(loci)}")
