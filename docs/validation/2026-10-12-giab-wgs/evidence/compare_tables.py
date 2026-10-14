#!/usr/bin/env python3
"""Compare the primary alignments of two aligners from bench/bam_summary.py tables.

Per read: identical (same contig, strand, start and end), same locus (same contig and strand,
overlapping), moved (different contig/strand or no overlap), newly unmapped / newly mapped,
and MAPQ changes (lower / higher; and reads crossing MAPQ 20 either way).

Usage: compare_tables.py REF/primary.tsv.gz TEST/primary.tsv.gz
"""
import gzip
import sys


def read(path):
    with gzip.open(path, "rt") as f:
        for line in f:
            name, flag, ctg, s, e, q = line.rstrip("\n").split("\t")
            yield name, int(flag), ctg, int(s), int(e), int(q)


ref = {}
for name, flag, ctg, s, e, q in read(sys.argv[1]):
    ref[name] = (flag, ctg, s, e, q)

c = dict(compared=0, identical=0, same_locus=0, moved=0, both_unmapped=0, newly_unmapped=0,
         newly_mapped=0, mapq_lower=0, mapq_higher=0, mapq_lost_20=0, mapq_gained_20=0, only_in_test=0)
for name, flag, ctg, s, e, q in read(sys.argv[2]):
    r = ref.pop(name, None)
    if r is None:
        c["only_in_test"] += 1
        continue
    c["compared"] += 1
    rflag, rctg, rs, re_, rq = r
    ru, tu = bool(rflag & 4), bool(flag & 4)
    if ru and tu:
        c["both_unmapped"] += 1
        continue
    if tu:
        c["newly_unmapped"] += 1
        continue
    if ru:
        c["newly_mapped"] += 1
        continue
    same_strand = (rflag & 16) == (flag & 16)
    if ctg == rctg and same_strand and s == rs and e == re_:
        c["identical"] += 1
    elif ctg == rctg and same_strand and s < re_ and rs < e:
        c["same_locus"] += 1
    else:
        c["moved"] += 1
    if q < rq:
        c["mapq_lower"] += 1
    elif q > rq:
        c["mapq_higher"] += 1
    if rq >= 20 > q:
        c["mapq_lost_20"] += 1
    elif q >= 20 > rq:
        c["mapq_gained_20"] += 1
c["only_in_ref"] = len(ref)

tot = max(c["compared"], 1)
for k, v in c.items():
    pct = f"\t{100 * v / tot:.3f}%" if k not in ("compared", "only_in_test", "only_in_ref") else ""
    print(f"{k}\t{v}{pct}")
