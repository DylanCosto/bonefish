#!/usr/bin/env python3
"""Write a synthetic large reference for limit tests: N random 10 kb contigs (decoy_*, which attract
no reads), then GRCh38 with chr1, chr2 and chr3 fused into one 689 Mbp contig (chr1_2_3, beyond the
2^29 positions a BAI index can address), then the other chromosomes unchanged. With 150,000
decoys the genome is 4.6 Gbp and chr19 onwards lie past 2^32 in the concatenated sequence.
Usage: large_ref.py GRCh38.fa N_DECOYS OUT.fa [random|N|pad]
With N decoys (all N bases) the coordinates shift the same way but the decoys add no seeds, so
alignments to GRCh38 must be unchanged; random decoys also add chance k-mer matches.
"pad" writes no decoys and leaves the contigs as they are, but puts 1.5 Gbp of N before chr1's
sequence: contig order (which minimap2-style aligners use to break ties) is unchanged, chr1 becomes a
1.75 Gbp contig, and every later chromosome lies 1.5 Gbp further along the concatenated reference.
"""
import random
import sys

src, n_decoy, out = sys.argv[1], int(sys.argv[2]), sys.argv[3]
kind = sys.argv[4] if len(sys.argv) > 4 else "random"
rng = random.Random(42)
table = bytes.maketrans(bytes(range(256)), b"ACGT" * 64)


def wrap(seq, fh):
    for i in range(0, len(seq), 80):
        fh.write(seq[i:i + 80]); fh.write(b"\n")


def read_fasta(path):
    name, chunks = None, []
    with open(path, "rb") as fh:
        for line in fh:
            if line.startswith(b">"):
                if name:
                    yield name, b"".join(chunks)
                name, chunks = line[1:].split()[0].decode(), []
            else:
                chunks.append(line.rstrip())
    if name:
        yield name, b"".join(chunks)


if kind == "pad":
    with open(out, "wb") as fh:
        for name, seq in read_fasta(src):
            fh.write(f">{name}\n".encode())
            wrap(b"N" * 1_500_000_000 + seq if name == "chr1" else seq, fh)
    sys.exit(0)

with open(out, "wb") as fh:
    for i in range(n_decoy):
        fh.write(f">decoy_{i:06d}\n".encode())
        wrap(rng.randbytes(10_000).translate(table) if kind == "random" else b"N" * 10_000, fh)
    fused = []
    for name, seq in read_fasta(src):
        if name in ("chr1", "chr2", "chr3"):
            fused.append(seq)
            if len(fused) == 3:
                fh.write(b">chr1_2_3\n"); wrap(b"".join(fused), fh); fused = None
            continue
        fh.write(f">{name}\n".encode()); wrap(seq, fh)
