#!/usr/bin/env python3
"""Compare spliced alignments of the same reads from two aligners, and against an annotation.

Per aligner (primary alignments): reads, mapped, spliced (any N in CIGAR), introns, and the share of
introns whose exact coordinates are annotated. Between aligners: reads whose primary intron chain
(ordered list of (contig, start, end) introns) is identical, among reads spliced by either.

Usage: intron_chains.py ANNOT.bed12 A.bam B.bam   (A = reference aligner, e.g. minimap2)
"""
import re
import subprocess
import sys

OP = re.compile(r"(\d+)([MIDNSHP=X])")


def annotated_introns(bed):
    s = set()
    for line in open(bed):
        f = line.rstrip("\n").split("\t")
        if len(f) < 12:
            continue
        c, start = f[0], int(f[1])
        sizes = [int(x) for x in f[10].rstrip(",").split(",")]
        starts = [int(x) for x in f[11].rstrip(",").split(",")]
        for i in range(len(sizes) - 1):
            s.add((c, start + starts[i] + sizes[i], start + starts[i + 1]))
    return s


def chains(bam):
    out = {}
    p = subprocess.Popen(["samtools", "view", "-F", "0x900", bam], stdout=subprocess.PIPE, text=True)
    for line in p.stdout:
        f = line.split("\t", 6)
        if int(f[1]) & 4:
            out[f[0]] = None
            continue
        pos = int(f[3]) - 1
        ch = []
        for n, op in OP.findall(f[5]):
            n = int(n)
            if op == "N":
                ch.append((f[2], pos, pos + n))
            if op in "MDN=X":
                pos += n
        out[f[0]] = tuple(ch)
    p.wait()
    return out


def stats(name, ch, annot):
    n = len(ch)
    mapped = sum(v is not None for v in ch.values())
    spliced = sum(bool(v) for v in ch.values())
    introns = [i for v in ch.values() if v for i in v]
    ann = sum(i in annot for i in introns)
    print(f"{name}: reads {n:,}  mapped {mapped:,} ({100 * mapped / max(n, 1):.2f}%)  "
          f"spliced {spliced:,} ({100 * spliced / max(mapped, 1):.2f}% of mapped)  "
          f"introns {len(introns):,}, annotated {100 * ann / max(len(introns), 1):.2f}%  "
          f"distinct introns {len(set(introns)):,}")


def main():
    annot = annotated_introns(sys.argv[1])
    a, b = chains(sys.argv[2]), chains(sys.argv[3])
    stats("A " + sys.argv[2], a, annot)
    stats("B " + sys.argv[3], b, annot)
    both = [q for q in a if q in b]
    spl = [q for q in both if a[q] or b[q]]
    same = sum(a[q] == b[q] for q in spl)
    only_a = sum(bool(a[q]) and not b[q] for q in spl)
    only_b = sum(bool(b[q]) and not a[q] for q in spl)
    print(f"reads spliced by either: {len(spl):,}; identical intron chain {same:,} "
          f"({100 * same / max(len(spl), 1):.2f}%); spliced only by A {only_a:,}, only by B {only_b:,}")
    # where chains differ, which one agrees better with the annotation
    a_better = b_better = 0
    for q in spl:
        if a[q] == b[q]:
            continue
        fa = sum(i in annot for i in (a[q] or ())) / max(len(a[q] or ()), 1)
        fb = sum(i in annot for i in (b[q] or ())) / max(len(b[q] or ()), 1)
        a_better += fa > fb
        b_better += fb > fa
    print(f"differing chains: higher annotated-intron share in A {a_better:,}, in B {b_better:,}")


if __name__ == "__main__":
    main()
