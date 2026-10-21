#!/usr/bin/env python3
"""Mapping statistics of one BAM for the metagenomics / duplex checks.

Prints reads, mapped share, MAPQ classes, supplementary records, aligned bases, gap-compressed
and BLAST-like identity of primaries (from NM), primaries carrying MM/ML tags, and, with a
circular-contig list, reads whose primary and a supplementary lie on the same circular contig
(origin-spanning reads). Writes per-contig primary counts (MAPQ >= 0 and >= 20) to COUNTS.tsv.

Usage: aln_stats.py ALN.bam COUNTS.tsv [circular.list]
"""
import collections
import re
import subprocess
import sys

OP = re.compile(r"(\d+)([MIDNSHP=X])")


def main():
    bam, counts_out = sys.argv[1], sys.argv[2]
    circ = set()
    if len(sys.argv) > 3:
        circ = {l.split()[0] for l in open(sys.argv[3]) if l.strip()}
    n = mapped = q0 = q20 = supp = mm = 0
    cols = nm_sum = gaps = 0
    aln_bases = 0
    per = collections.defaultdict(lambda: [0, 0])
    prim_ctg = {}
    supp_same_circ = set()
    pending = []
    p = subprocess.Popen(["samtools", "view", "-F", "0x100", bam], stdout=subprocess.PIPE, text=True)
    for line in p.stdout:
        f = line.rstrip("\n").split("\t")
        flag = int(f[1])
        if flag & 0x800:
            supp += 1
            pending.append((f[0], f[2]))
            continue
        n += 1
        if flag & 0x4:
            continue
        mapped += 1
        mapq = int(f[4])
        q0 += mapq == 0
        q20 += mapq >= 20
        per[f[2]][0] += 1
        per[f[2]][1] += mapq >= 20
        prim_ctg[f[0]] = f[2]
        m = ins = dele = nins = ndel = 0
        for k, op in OP.findall(f[5]):
            k = int(k)
            if op in "M=X":
                m += k
            elif op == "I":
                ins += k
                nins += 1
            elif op == "D":
                dele += k
                ndel += 1
        aln_bases += m
        nm = next((int(t[5:]) for t in f[11:] if t.startswith("NM:i:")), None)
        if nm is not None:
            cols += m + ins + dele
            nm_sum += nm
            # gap-compressed: each indel counts once
            gaps += nm - ins - dele + nins + ndel
        mm += any(t.startswith("MM:Z:") or t.startswith("Mm:Z:") for t in f[11:])
    p.wait()
    for q, c in pending:
        if c in circ and prim_ctg.get(q) == c:
            supp_same_circ.add(q)
    print(f"reads {n:,}  mapped {mapped:,} ({100 * mapped / max(n, 1):.2f}%)  "
          f"MAPQ0 {q0:,}  MAPQ>=20 {q20:,} ({100 * q20 / max(n, 1):.2f}%)")
    print(f"supplementary {supp:,}  aligned bases (primary) {aln_bases / 1e9:.3f} Gb")
    if cols:
        gc_cols = cols - (nm_sum - gaps)
        print(f"identity BLAST-like {1 - nm_sum / cols:.5f}  gap-compressed {1 - gaps / max(gc_cols, 1):.5f}")
    print(f"primaries with MM tag {mm:,} ({100 * mm / max(mapped, 1):.2f}% of mapped)")
    if circ:
        print(f"origin-spanning reads on circular contigs (primary + supplementary on same contig): {len(supp_same_circ):,}")
    with open(counts_out, "w") as o:
        for c in sorted(per):
            o.write(f"{c}\t{per[c][0]}\t{per[c][1]}\n")


if __name__ == "__main__":
    main()
