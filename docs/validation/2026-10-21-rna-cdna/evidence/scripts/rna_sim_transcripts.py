#!/usr/bin/env python3
"""Transcript sequences for simulated spliced reads with known intron chains.

Picks N random multi-exon transcripts (protein-coding and lncRNA, length >= 300 bp) from a GENCODE
BED12 and writes their mRNA sequences (reverse-complemented for minus-strand transcripts, so every
sequence is 5'->3') plus a truth table of genomic exon blocks. Badread then samples reads from these
sequences; rna_sim_truth.py converts read coordinates back to genomic introns.

Usage: rna_sim_transcripts.py GENOME.fa(+.fai) ANNOT.bed12 N SEED OUT_PREFIX
"""
import mmap
import random
import sys

COMP = bytes.maketrans(b"ACGTNacgtn", b"TGCANtgcan")


def main():
    fa, bed, n, seed, out = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4]), sys.argv[5]
    fai = {}
    for line in open(fa + ".fai"):
        c, ln, off, lb, lby = line.split()[:5]
        fai[c] = (int(ln), int(off), int(lb), int(lby))
    tx = []
    for line in open(bed):
        f = line.rstrip("\n").split("\t")
        if len(f) < 12 or f[0] not in fai or int(f[9]) < 2:
            continue
        sizes = [int(x) for x in f[10].rstrip(",").split(",")]
        if sum(sizes) < 300:
            continue
        tx.append(f)
    rng = random.Random(seed)
    pick = rng.sample(tx, min(n, len(tx)))
    fh = open(fa, "rb")
    mm = mmap.mmap(fh.fileno(), 0, access=mmap.ACCESS_READ)

    def seq(c, s, e):
        ln, off, lb, lby = fai[c]
        a = off + (s // lb) * lby + s % lb
        b = off + (e // lb) * lby + e % lb
        return mm[a:b].replace(b"\n", b"").upper()

    with open(out + ".fa", "w") as fo, open(out + ".truth.tsv", "w") as ft:
        for f in pick:
            c, start, name, strand = f[0], int(f[1]), f[3], f[5]
            sizes = [int(x) for x in f[10].rstrip(",").split(",")]
            starts = [int(x) for x in f[11].rstrip(",").split(",")]
            blocks = [(start + s, start + s + z) for s, z in zip(starts, sizes)]
            s = b"".join(seq(c, a, b) for a, b in blocks)
            if strand == "-":
                s = s.translate(COMP)[::-1]
            fo.write(f">{name}\n{s.decode()}\n")
            ft.write(f"{name}\t{c}\t{strand}\t" + ",".join(f"{a}-{b}" for a, b in blocks) + "\n")


if __name__ == "__main__":
    main()
