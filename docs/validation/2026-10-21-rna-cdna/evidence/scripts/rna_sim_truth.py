#!/usr/bin/env python3
"""Score spliced alignments of Badread reads simulated from transcripts (rna_sim_transcripts.py).

Badread read names carry `transcript,strand,start-end` (transcript coordinates of the read). The true
genomic intron chain of a read is the set of introns between exon blocks that the read's transcript
interval spans (a read must cover >= 10 bp on both sides of an intron for it to count). Per aligner:
reads mapped, reads with the exact true intron chain, intron precision / recall, wrong-locus primaries.

Usage: rna_sim_truth.py TRUTH.tsv READS.fq.gz A.bam [B.bam ...]
"""
import gzip
import re
import subprocess
import sys

OP = re.compile(r"(\d+)([MIDNSHP=X])")
MIN_FLANK = 10


def load_truth(path):
    t = {}
    for line in open(path):
        name, c, strand, bl = line.rstrip("\n").split("\t")
        blocks = [tuple(map(int, b.split("-"))) for b in bl.split(",")]
        t[name] = (c, strand, blocks)
    return t


def true_chain(t, desc):
    # Badread description: "TRANSCRIPT,+strand,START-END length=... error-free_length=... read_identity=..."
    m = re.match(r"(\S+),([+-])strand,(\d+)-(\d+)", desc)
    if not m:
        return None
    name, s, e = m.group(1), int(m.group(3)), int(m.group(4))
    c, strand, blocks = t[name]
    L = sum(b - a for a, b in blocks)
    # transcript coordinates are 5'->3'; convert to genomic-forward exonic offsets
    if strand == "-":
        s, e = L - e, L - s
    introns = []
    off = 0
    for i, (a, b) in enumerate(blocks):
        off += b - a
        if i + 1 < len(blocks) and s + MIN_FLANK <= off <= e - MIN_FLANK:
            introns.append((c, b, blocks[i + 1][0]))
    # genomic span of the read
    gs = ge = None
    off = 0
    for a, b in blocks:
        if gs is None and s < off + (b - a):
            gs = a + (s - off)
        if e <= off + (b - a):
            ge = a + (e - off)
            break
        off += b - a
    return c, gs, ge, tuple(introns)


def chains(bam):
    out = {}
    p = subprocess.Popen(["samtools", "view", "-F", "0x900", bam], stdout=subprocess.PIPE, text=True)
    for line in p.stdout:
        f = line.split("\t", 6)
        if int(f[1]) & 4:
            out[f[0]] = None
            continue
        pos = int(f[3]) - 1
        start = pos
        ch = []
        for n, op in OP.findall(f[5]):
            n = int(n)
            if op == "N":
                ch.append((f[2], pos, pos + n))
            if op in "MDN=X":
                pos += n
        out[f[0]] = (f[2], start, pos, tuple(ch))
    p.wait()
    return out


def main():
    t = load_truth(sys.argv[1])
    truth = {}
    with gzip.open(sys.argv[2], "rt") as fq:
        for i, line in enumerate(fq):
            if i % 4 == 0:
                rid, _, desc = line[1:].rstrip("\n").partition(" ")
                tc = true_chain(t, desc)
                if tc:
                    truth[rid] = tc
    n = len(truth)
    spliced_truth = sum(bool(v[3]) for v in truth.values())
    print(f"simulated reads {n:,}, with >= 1 true intron {spliced_truth:,}")
    for bam in sys.argv[3:]:
        a = chains(bam)
        mapped = exact = wrong = tp = fp = fn = 0
        for rid, (c, gs, ge, tch) in truth.items():
            r = a.get(rid)
            if r is None:
                fn += len(tch)
                continue
            mapped += 1
            if r[0] != c or r[2] < gs or r[1] > ge:
                wrong += 1
            exact += r[3] == tch
            ts, rs = set(tch), set(r[3])
            tp += len(ts & rs)
            fp += len(rs - ts)
            fn += len(ts - rs)
        print(f"{bam}: mapped {100 * mapped / n:.2f}%  exact intron chain {exact:,} ({100 * exact / n:.2f}%)  "
              f"wrong locus {wrong:,}  intron precision {100 * tp / max(tp + fp, 1):.2f}%  "
              f"recall {100 * tp / max(tp + fn, 1):.2f}%")


if __name__ == "__main__":
    main()
