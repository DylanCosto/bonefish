#!/usr/bin/env python3
"""Compare per-contig and per-bin read counts of two aligners (metagenomics check).

Inputs are aln_stats.py COUNTS.tsv files (contig, primaries, primaries with MAPQ >= 20) and a
contig-to-bin table (DAS Tool). Prints Pearson correlation of counts and of log counts, the
largest relative differences among contigs/bins with >= 100 reads, and the share of reads whose
contig or bin assignment would differ in total.

Usage: abundance_compare.py A.tsv B.tsv contig2bin.tsv
"""
import math
import sys


def load(path):
    d = {}
    for line in open(path):
        c, n, q = line.split()
        d[c] = (int(n), int(q))
    return d


def pearson(x, y):
    n = len(x)
    mx, my = sum(x) / n, sum(y) / n
    sxy = sum((a - mx) * (b - my) for a, b in zip(x, y))
    sx = math.sqrt(sum((a - mx) ** 2 for a in x))
    sy = math.sqrt(sum((b - my) ** 2 for b in y))
    return sxy / (sx * sy) if sx and sy else float("nan")


def report(name, a, b, col):
    keys = sorted(set(a) | set(b))
    x = [a.get(k, (0, 0))[col] for k in keys]
    y = [b.get(k, (0, 0))[col] for k in keys]
    tot = sum(x)
    diff = sum(abs(p - q) for p, q in zip(x, y)) / 2
    lx = [math.log10(v + 1) for v in x]
    ly = [math.log10(v + 1) for v in y]
    print(f"{name}: {len(keys):,} with reads; total {tot:,} vs {sum(y):,}; "
          f"Pearson {pearson(x, y):.6f}, on log10(n+1) {pearson(lx, ly):.6f}; "
          f"net reassigned {int(diff):,} ({100 * diff / max(tot, 1):.3f}%)")
    big = [(abs(q - p) / max(p, 1), k, p, q) for k, p, q in zip(keys, x, y) if max(p, q) >= 100]
    big.sort(reverse=True)
    for r, k, p, q in big[:5]:
        print(f"    {k}: {p:,} vs {q:,} ({100 * (q - p) / max(p, 1):+.1f}%)")


def main():
    a, b = load(sys.argv[1]), load(sys.argv[2])
    c2b = {}
    for line in open(sys.argv[3]):
        f = line.split()
        if len(f) >= 2:
            c2b[f[0]] = f[1]
    for col, label in ((0, "all primaries"), (1, "MAPQ >= 20")):
        print(f"== {label}")
        report("contigs", a, b, col)
        ba, bb = {}, {}
        for src, dst in ((a, ba), (b, bb)):
            for c, v in src.items():
                k = c2b.get(c, "unbinned")
                n, q = dst.get(k, (0, 0))
                dst[k] = (n + v[0], q + v[1])
        report("bins", ba, bb, col)


if __name__ == "__main__":
    main()
