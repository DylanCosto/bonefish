#!/usr/bin/env python3
"""One-pass summary of a whole-genome BAM, from `samtools view BAM | cut -f1-6,12-` on stdin
(SEQ/QUAL dropped so the pass is cheap).

Writes OUT_DIR/bam_summary.tsv (record counts, primary query bases and identity, MAPQ
distribution) and OUT_DIR/primary.tsv.gz (one row per read: name, flag, contig, start, end,
mapq) for comparing aligners with bench/compare_tables.py.

Identity is BLAST-like over primary alignments: 1 - NM / (M + I + D columns).

Usage: samtools view -@4 aln.bam | cut -f1-6,12- | bam_summary.py OUT_DIR
"""
import gzip
import re
import sys

out_dir = sys.argv[1]
CIG = re.compile(r"(\d+)([MIDNSHP=X])")
MAPQ_BINS = [(0, 0), (1, 9), (10, 19), (20, 29), (30, 59), (60, 60)]

n = {"records": 0, "primary_mapped": 0, "unmapped": 0, "secondary": 0, "supplementary": 0}
q_primary = 0          # query bases of primary records (mapped + unmapped via CIGAR/none)
aligned_q = 0          # query bases inside primary alignments (M/I/=/X)
cols = 0               # M + I + D alignment columns, primary
nm = 0
supp_aligned_q = 0
mapq_hist = [0] * len(MAPQ_BINS)

with gzip.open(f"{out_dir}/primary.tsv.gz", "wt", compresslevel=1) as tab:
    for line in sys.stdin:
        f = line.rstrip("\n").split("\t")
        flag = int(f[1])
        n["records"] += 1
        if flag & 0x4:
            n["unmapped"] += 1
            tab.write(f"{f[0]}\t{flag}\t*\t0\t0\t0\n")
            continue
        if flag & 0x100:
            n["secondary"] += 1
            continue
        m = i = d = s = 0
        for ln, op in CIG.findall(f[5]):
            ln = int(ln)
            if op in "M=X":
                m += ln
            elif op == "I":
                i += ln
            elif op in "DN":
                d += ln
            elif op == "S":
                s += ln
        if flag & 0x800:
            n["supplementary"] += 1
            supp_aligned_q += m + i
            continue
        n["primary_mapped"] += 1
        q_primary += m + i + s
        aligned_q += m + i
        cols += m + i + d
        for t in f[6:]:
            if t.startswith("NM:i:"):
                nm += int(t[5:])
                break
        mapq = int(f[4])
        for b, (lo, hi) in enumerate(MAPQ_BINS):
            if lo <= mapq <= hi:
                mapq_hist[b] += 1
                break
        start = int(f[3]) - 1
        tab.write(f"{f[0]}\t{flag}\t{f[2]}\t{start}\t{start + m + d}\t{mapq}\n")

with open(f"{out_dir}/bam_summary.tsv", "w") as o:
    for k, v in n.items():
        o.write(f"{k}\t{v}\n")
    reads = n["primary_mapped"] + n["unmapped"]
    o.write(f"reads\t{reads}\n")
    o.write(f"pct_reads_mapped\t{100 * n['primary_mapped'] / max(reads, 1):.3f}\n")
    o.write(f"primary_query_bases\t{q_primary}\n")
    o.write(f"primary_aligned_bases\t{aligned_q}\n")
    o.write(f"pct_primary_bases_aligned\t{100 * aligned_q / max(q_primary, 1):.3f}\n")
    o.write(f"supplementary_aligned_bases\t{supp_aligned_q}\n")
    o.write(f"primary_identity\t{1 - nm / max(cols, 1):.5f}\n")
    for (lo, hi), c in zip(MAPQ_BINS, mapq_hist):
        o.write(f"mapq_{lo}-{hi}\t{c}\n")
    o.write(f"pct_primary_mapq_ge_20\t{100 * sum(mapq_hist[3:]) / max(n['primary_mapped'], 1):.3f}\n")
