#!/usr/bin/env python3
"""Compare alignments of two BAMs of the same reads (Phase 3 output-change audit).

Primary records, over reads present in the reference BAM:
  identical              same contig, start, strand, CIGAR
  same locus             same contig and strand, reference intervals overlap
  moved                  different contig/strand, or no overlap
  newly unmapped/mapped
Also totals of aligned query bases (M/=/X) for primary and supplementary records and of
soft-clipped bases on primaries, to see whether sequence cut from one alignment is recovered
as a supplementary alignment.

Usage: compare_primary.py REF.bam TEST.bam
"""
import re
import subprocess
import sys

OP = re.compile(r"(\d+)([MIDNSHP=X])")


def cigar_stats(cig):
    aln = clip = rlen = 0
    for n, op in OP.findall(cig):
        n = int(n)
        if op in "M=X":
            aln += n
            rlen += n
        elif op in "DN":
            rlen += n
        elif op == "S":
            clip += n
    return aln, clip, rlen


def load(bam):
    prim, tot = {}, {"prim_aln": 0, "prim_clip": 0, "supp_aln": 0, "supp_n": 0}
    p = subprocess.Popen(["samtools", "view", "-F", "0x100", bam], stdout=subprocess.PIPE, text=True)
    for line in p.stdout:
        f = line.split("\t", 6)
        flag = int(f[1])
        if flag & 0x4:
            prim[f[0]] = (False, None, 0, 0, False, 0, "*")
            continue
        aln, clip, rlen = cigar_stats(f[5])
        if flag & 0x800:
            tot["supp_aln"] += aln
            tot["supp_n"] += 1
            continue
        tot["prim_aln"] += aln
        tot["prim_clip"] += clip
        start = int(f[3])
        prim[f[0]] = (True, f[2], start, start + rlen, bool(flag & 0x10), int(f[4]), f[5])
    p.wait()
    return prim, tot


def main():
    ref, rt = load(sys.argv[1])
    test, tt = load(sys.argv[2])
    n = len(ref)
    same = locus = moved = lost = gained = mapq_down = mapq_up = 0
    for q, r in ref.items():
        t = test.get(q)
        if t is None:
            continue
        if r[0] and not t[0]:
            lost += 1
        elif not r[0] and t[0]:
            gained += 1
        elif not r[0]:
            same += 1
        elif (r[1], r[2], r[4], r[6]) == (t[1], t[2], t[4], t[6]):
            same += 1
        elif r[1] == t[1] and r[4] == t[4] and r[2] < t[3] and t[2] < r[3]:
            locus += 1
        else:
            moved += 1
        if r[0] and t[0]:
            mapq_down += t[5] < r[5]
            mapq_up += t[5] > r[5]
    pct = lambda x: f"{x:,} ({100 * x / max(n, 1):.3f}%)"
    gb = lambda x: f"{x / 1e9:.3f} Gb"
    print(f"primary records compared: {n:,}")
    print(f"  identical:                         {pct(same)}")
    print(f"  same locus (overlapping interval): {pct(locus)}")
    print(f"  moved (contig/strand/no overlap):  {pct(moved)}")
    print(f"  newly unmapped / newly mapped:     {lost:,} / {gained:,}")
    print(f"  MAPQ lower / higher:               {mapq_down:,} / {mapq_up:,}")
    print(f"aligned bases, primary:        {gb(rt['prim_aln'])} -> {gb(tt['prim_aln'])}")
    print(f"aligned bases, supplementary:  {gb(rt['supp_aln'])} -> {gb(tt['supp_aln'])} "
          f"({rt['supp_n']:,} -> {tt['supp_n']:,} records)")
    print(f"aligned bases, total:          {gb(rt['prim_aln'] + rt['supp_aln'])} -> {gb(tt['prim_aln'] + tt['supp_aln'])}")
    print(f"soft-clipped bases, primary:   {gb(rt['prim_clip'])} -> {gb(tt['prim_clip'])}")


if __name__ == "__main__":
    main()
