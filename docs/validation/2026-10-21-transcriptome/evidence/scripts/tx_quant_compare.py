#!/usr/bin/env python3
"""Compare oarfish quantifications from two aligners' transcriptome BAMs.

Reads `<prefix>.quant` (tname, len, num_reads) of each, and prints Pearson correlation of counts and
of log2(count + 1), Spearman correlation, and the share of reads assigned differently (half the sum
of absolute count differences / total), at transcript and gene level (gene = 2nd |-field of the
GENCODE transcript name). Only transcripts with at least one read in either run are compared.

Usage: tx_quant_compare.py A.quant B.quant
"""
import math
import sys


def load(path):
    d = {}
    with open(path) as f:
        header = f.readline().rstrip("\n").split("\t")
        ci = header.index("num_reads") if "num_reads" in header else 2
        for line in f:
            x = line.rstrip("\n").split("\t")
            d[x[0]] = float(x[ci])
    return d


def pearson(x, y):
    n = len(x)
    mx, my = sum(x) / n, sum(y) / n
    sxy = sum((a - mx) * (b - my) for a, b in zip(x, y))
    sx = math.sqrt(sum((a - mx) ** 2 for a in x))
    sy = math.sqrt(sum((b - my) ** 2 for b in y))
    return sxy / (sx * sy) if sx and sy else float("nan")


def ranks(v):
    order = sorted(range(len(v)), key=lambda i: v[i])
    r = [0.0] * len(v)
    i = 0
    while i < len(v):
        j = i
        while j + 1 < len(v) and v[order[j + 1]] == v[order[i]]:
            j += 1
        for k in range(i, j + 1):
            r[order[k]] = (i + j) / 2 + 1
        i = j + 1
    return r


def report(name, a, b):
    keys = [k for k in set(a) | set(b) if a.get(k, 0) > 0 or b.get(k, 0) > 0]
    x = [a.get(k, 0.0) for k in keys]
    y = [b.get(k, 0.0) for k in keys]
    tot = sum(x)
    diff = sum(abs(p - q) for p, q in zip(x, y)) / 2
    lx = [math.log2(v + 1) for v in x]
    ly = [math.log2(v + 1) for v in y]
    print(f"{name}: {len(keys):,} with reads; reads {tot:,.0f} vs {sum(y):,.0f}; Pearson {pearson(x, y):.5f}, "
          f"log2 {pearson(lx, ly):.5f}, Spearman {pearson(ranks(x), ranks(y)):.5f}; "
          f"reassigned {diff:,.0f} ({100 * diff / max(tot, 1):.2f}%)")


def gene_level(d):
    g = {}
    for k, v in d.items():
        f = k.split("|")
        gid = f[1] if len(f) > 1 else k
        g[gid] = g.get(gid, 0.0) + v
    return g


def main():
    a, b = load(sys.argv[1]), load(sys.argv[2])
    report("transcripts", a, b)
    report("genes", gene_level(a), gene_level(b))


if __name__ == "__main__":
    main()
