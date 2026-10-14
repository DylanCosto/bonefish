#!/usr/bin/env python3
"""Per-site comparison of two modkit bedMethyl files (e.g. CpG pileups from two aligners' BAMs).

For each modification code: sites with valid coverage >= MIN_COV in both files, Pearson r and
mean/median absolute difference of the percent modified, the share of sites differing by > 10 and
> 25 points, coverage agreement, sites covered (>= MIN_COV) in only one file, the overall level
(sum modified / sum valid), and the correlation of 1 kb bins (coverage-weighted).

Usage: modkit_compare.py A.bed B.bed [NAME_A NAME_B] [MIN_COV]
"""
import sys
import numpy as np

a_path, b_path = sys.argv[1], sys.argv[2]
na, nb = (sys.argv[3], sys.argv[4]) if len(sys.argv) > 4 else ("A", "B")
min_cov = int(sys.argv[5]) if len(sys.argv) > 5 else 10


def load(path):
    d = {}
    for line in open(path):
        f = line.split()
        key = (f[0], int(f[1]), f[3])
        d[key] = (int(f[9]), float(f[10]), int(f[11]))  # valid coverage, percent modified, n modified
    return d


A, B = load(a_path), load(b_path)
codes = sorted({k[2] for k in A} | {k[2] for k in B})
print(f"{na} vs {nb}, sites with valid coverage >= {min_cov}")
print("code\tsites_A\tsites_B\tcommon\tonly_A\tonly_B\tpearson_r\tmean_abs_diff\tmedian_abs_diff\tpct_diff_gt10\tpct_diff_gt25\t"
      "level_A\tlevel_B\tmean_cov_A\tmean_cov_B\tcov_r\tbin1k_r")
for c in codes:
    ka = {k for k, v in A.items() if k[2] == c and v[0] >= min_cov}
    kb = {k for k, v in B.items() if k[2] == c and v[0] >= min_cov}
    com = sorted(ka & kb)
    pa = np.array([A[k][1] for k in com]); pb = np.array([B[k][1] for k in com])
    ca = np.array([A[k][0] for k in com]); cb = np.array([B[k][0] for k in com])
    d = np.abs(pa - pb)
    lev = lambda D: sum(v[2] for k, v in D.items() if k[2] == c) / max(1, sum(v[0] for k, v in D.items() if k[2] == c))
    # 1 kb bins over all sites of this code (coverage-weighted level per bin), bins with >= 50 valid calls in both
    def bins(D):
        m, n = {}, {}
        for k, v in D.items():
            if k[2] != c: continue
            b = (k[0], k[1] // 1000)
            m[b] = m.get(b, 0) + v[2]; n[b] = n.get(b, 0) + v[0]
        return m, n
    ma, nA = bins(A); mb, nB = bins(B)
    bk = [b for b in nA if nA[b] >= 50 and nB.get(b, 0) >= 50]
    ba = np.array([ma[b] / nA[b] for b in bk]); bb = np.array([mb[b] / nB[b] for b in bk])
    r = np.corrcoef(pa, pb)[0, 1] if len(com) > 2 else float("nan")
    rc = np.corrcoef(ca, cb)[0, 1] if len(com) > 2 else float("nan")
    rb = np.corrcoef(ba, bb)[0, 1] if len(bk) > 2 else float("nan")
    print(f"{c}\t{len(ka)}\t{len(kb)}\t{len(com)}\t{len(ka - kb)}\t{len(kb - ka)}\t{r:.4f}\t{d.mean():.2f}\t{np.median(d):.2f}\t"
          f"{100 * (d > 10).mean():.2f}\t{100 * (d > 25).mean():.2f}\t{100 * lev(A):.2f}\t{100 * lev(B):.2f}\t"
          f"{ca.mean():.1f}\t{cb.mean():.1f}\t{rc:.4f}\t{rb:.4f}")
