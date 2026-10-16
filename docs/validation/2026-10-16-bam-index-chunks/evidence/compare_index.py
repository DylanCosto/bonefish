#!/usr/bin/env python3
"""Compare two BAM indexes (both BAI, or both CSI) by content: bins, chunks and linear index/loff
per reference (the order bins are stored in is ignored). Usage: compare_index.py A.bai B.bai"""
import gzip, struct, sys

def load(p):
    d = open(p, 'rb').read()
    csi = d[:2] == b'\x1f\x8b'
    if csi: d = gzip.decompress(d)
    o = 4
    if csi:
        ms, depth, laux = struct.unpack_from('<iii', d, o); o += 12 + laux
    (n_ref,) = struct.unpack_from('<i', d, o); o += 4
    refs = []
    for _ in range(n_ref):
        (n_bin,) = struct.unpack_from('<i', d, o); o += 4
        bins = {}
        for _ in range(n_bin):
            if csi: b, loff, nc = struct.unpack_from('<IQi', d, o); o += 16
            else: (b, nc), loff = struct.unpack_from('<Ii', d, o), None; o += 8
            bins[b] = (loff, [struct.unpack_from('<QQ', d, o + 16 * i) for i in range(nc)]); o += 16 * nc
        lin = []
        if not csi:
            (n_intv,) = struct.unpack_from('<i', d, o); o += 4
            lin = list(struct.unpack_from('<%dQ' % n_intv, d, o)); o += 8 * n_intv
        refs.append((bins, lin))
    return ('CSI' if csi else 'BAI'), refs, d[o:]

ka, A, ta = load(sys.argv[1]); kb, B, tb = load(sys.argv[2])
chunks = lambda R: sum(len(v[1]) for bins, _ in R for k, v in bins.items() if k not in (37450,))
print(f'{ka} {sys.argv[1]}: {chunks(A)} chunks; {kb} {sys.argv[2]}: {chunks(B)} chunks')
print('identical content:', A == B and ta == tb)
