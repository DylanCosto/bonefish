#!/usr/bin/env python3
"""Mapping accuracy on Badread-simulated reads (true origin in the read name).

Reads with a single source segment ("chr20,+strand,START-END") are scored; chimeric reads (several
segments) are skipped. A primary alignment is correct if it is on the true contig and strand and
overlaps the true interval. Badread reports "-strand" coordinates on the reverse-complemented
contig, so they are converted with the contig length from @SQ. Reports, per MAPQ threshold, the
number of reads mapped at or above it and the fraction of those that are wrong (like minimap2's
paftools mapeval), the aligned fraction of read bases, and how many junk/random reads (no true
origin) were mapped, which should be none at high MAPQ.

Usage: bonefish -a ... reads.fq.gz | sim_eval.py LABEL reads.fq.gz [CURVE.tsv]
       (SAM on stdin; truths are read from the FASTQ comments, which the aligner does not copy
       without -y; CURVE.tsv gets mapped/wrong counts at every MAPQ for plotting)
"""
import re
import sys

OP = re.compile(r"(\d+)([MIDNSHP=X])")
TRUTH = re.compile(r"^([^,\s]+),([+-])strand,(\d+)-(\d+)$")


def main():
    label = sys.argv[1] if len(sys.argv) > 1 else "run"
    rows = []          # (mapq, correct) for mapped primaries
    n_scored = n_unmapped = n_noorigin = 0
    noorigin_mapq = []
    aln_bases = read_bases = 0
    lengths = {}
    for line in sys.stdin:
        if line.startswith("@"):
            if line.startswith("@SQ"):
                tags = dict(t.split(":", 1) for t in line.rstrip("\n").split("\t")[1:] if ":" in t)
                lengths[tags["SN"]] = int(tags["LN"])
            continue
        f = line.split("\t", 12)
        flag = int(f[1])
        if flag & 0x900:
            continue
        t = TRUTHS.get(f[0])
        if t is None:
            if f[0] in NOORIGIN:
                n_noorigin += 1
                if not flag & 0x4:
                    noorigin_mapq.append(int(f[4]))
            continue
        n_scored += 1
        read_bases += t[4]
        if flag & 0x4:
            n_unmapped += 1
            continue
        mapq = int(f[4])
        pos = int(f[3])
        rlen = 0
        for n, op in OP.findall(f[5]):
            if op in "MDN=X":
                rlen += int(n)
            if op in "M=X":
                aln_bases += int(n)
        strand = "-" if flag & 0x10 else "+"
        ts, te = t[2], t[3]
        if t[1] == "-":
            L = lengths[t[0]]
            ts, te = L - t[3], L - t[2]
        ok = f[2] == t[0] and strand == t[1] and pos <= te and pos + rlen >= ts
        rows.append((mapq, ok))
    print(f"## {label}: {n_scored:,} single-origin reads scored, {n_unmapped:,} unmapped, "
          f"primary aligned fraction {aln_bases / max(read_bases, 1):.4f}")
    print("| MAPQ >= | mapped | wrong | error rate |")
    print("|---:|---:|---:|---:|")
    for thr in (0, 1, 5, 10, 20, 30, 40, 50, 60):
        sel = [ok for q, ok in rows if q >= thr]
        wrong = sum(1 for ok in sel if not ok)
        print(f"| {thr} | {len(sel):,} | {wrong:,} | {wrong / max(len(sel), 1):.5f} |")
    print(f"junk/random reads: {n_noorigin:,}, mapped {len(noorigin_mapq):,}, "
          f"at MAPQ >= 10: {sum(q >= 10 for q in noorigin_mapq):,}, "
          f"at MAPQ >= 30: {sum(q >= 30 for q in noorigin_mapq):,}")
    if len(sys.argv) > 3:
        with open(sys.argv[3], "w") as out:
            out.write("mapq\tmapped_at_or_above\twrong_at_or_above\tscored\n")
            for thr in range(max((q for q, _ in rows), default=0), -1, -1):
                sel = [ok for q, ok in rows if q >= thr]
                out.write(f"{thr}\t{len(sel)}\t{sum(1 for ok in sel if not ok)}\t{n_scored}\n")


def load_truths(fastq):
    import gzip
    truths = {}
    with gzip.open(fastq, "rt") as fh:
        for i, line in enumerate(fh):
            if i % 4:
                continue
            parts = line[1:].split()
            if len(parts) < 3:
                continue
            if parts[1] in ("junk_seq", "random_seq"):
                NOORIGIN.add(parts[0])
                continue
            m = TRUTH.match(parts[1])
            if not m or "chimera" in parts:
                continue   # chimeric
            length = int(parts[2].split("=")[1]) if parts[2].startswith("length=") else 0
            truths[parts[0]] = (m.group(1), m.group(2), int(m.group(3)), int(m.group(4)), length)
    return truths


NOORIGIN = set()

if __name__ == "__main__":
    TRUTHS = load_truths(sys.argv[2])
    main()
