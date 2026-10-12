#!/usr/bin/env python3
"""Generate the small synthetic ONT-like dataset used by bonefish/tests/regression.rs.

Deterministic (fixed seed). Writes bonefish/tests/data/{ref.fa.gz,reads.fq.gz}. The reference
has three contigs with interspersed repeats, tandem repeats and a segmental duplication; the
reads (2-14% error, ONT-style homopolymer indels, both strands) include chimeras, reads spanning
a large deletion / insertion / inversion relative to the reference, and unmappable junk, so
the test exercises extension, gap filling, the RMQ rescue, supplementary/secondary output and
the unmapped path.

Expected outputs are produced separately (see the test's doc comment).
"""
import gzip
import os
import random

R = random.Random(20261006)
OUT = os.path.join(os.path.dirname(__file__), "..", "bonefish", "tests", "data")
COMP = str.maketrans("ACGT", "TGCA")


def rand_seq(n):
    return "".join(R.choice("ACGT") for _ in range(n))


def revcomp(s):
    return s.translate(COMP)[::-1]


def mutate(s, rate):
    """Copy with substitutions/indels at `rate` (for diverged repeat copies)."""
    out = []
    for c in s:
        x = R.random()
        if x < rate * 0.6:
            out.append(R.choice("ACGT"))
        elif x < rate * 0.8:
            continue
        elif x < rate:
            out.append(c + R.choice("ACGT"))
        else:
            out.append(c)
    return "".join(out)


def make_contig(n, repeat_pool):
    s = list(rand_seq(n))
    # interspersed repeats (diverged copies of a few repeat families)
    for _ in range(n // 20000):
        fam = R.choice(repeat_pool)
        pos = R.randrange(0, n - len(fam))
        s[pos:pos + len(fam)] = mutate(fam, 0.08)
    # tandem repeats
    for _ in range(n // 50000):
        unit = rand_seq(R.choice([2, 3, 6, 15, 40]))
        pos = R.randrange(0, n - 3000)
        tr = (unit * (R.randrange(400, 2500) // len(unit) + 1))
        s[pos:pos + len(tr)] = mutate(tr, 0.03)
    return "".join(s)


repeat_pool = [rand_seq(R.choice([300, 1000, 6000])) for _ in range(4)]
ref = {
    "ctgA": make_contig(400_000, repeat_pool),
    "ctgB": make_contig(200_000, repeat_pool),
    "ctgC": make_contig(60_000, repeat_pool),
}
# segmental duplication: 15 kb of ctgA copied (98% identity) into ctgB
a = ref["ctgA"]
dup = mutate(a[100_000:115_000], 0.02)
b = ref["ctgB"]
ref["ctgB"] = b[:120_000] + dup + b[120_000 + len(dup):]


def ont_errors(s, rate=0.04):
    out = []
    i = 0
    while i < len(s):
        c = s[i]
        # homopolymer run: ONT under/over-calls their length
        j = i
        while j < len(s) and s[j] == c:
            j += 1
        run = j - i
        if run >= 3 and R.random() < 0.25:
            run += R.choice([-1, -1, 1])
        for _ in range(run):
            x = R.random()
            if x < rate * 0.4:
                out.append(R.choice("ACGT"))
            elif x < rate * 0.7:
                pass
            elif x < rate:
                out.append(c)
                out.append(R.choice("ACGT"))
            else:
                out.append(c)
        i = j
    return "".join(out)


def qual(n):
    return "".join(chr(33 + max(3, min(50, int(R.gauss(22, 6))))) for _ in range(n))


def sample(length):
    name = R.choices(list(ref), weights=[len(v) for v in ref.values()])[0]
    seq = ref[name]
    length = min(length, len(seq) - 1)
    st = R.randrange(0, len(seq) - length)
    return seq[st:st + length]


reads = []
for i in range(160):
    L = int(min(40_000, max(800, R.lognormvariate(8.4, 0.7))))
    kind = R.random()
    if kind < 0.06:  # chimera
        s = sample(L // 2) + sample(L // 2)
    elif kind < 0.12:  # large deletion relative to the reference
        s = sample(L + 3000)
        cut = len(s) // 2
        s = s[:cut] + s[cut + R.choice([600, 3000, 8000]) % (len(s) - cut - 200):]
    elif kind < 0.17:  # large insertion
        s = sample(L)
        cut = len(s) // 2
        s = s[:cut] + rand_seq(R.choice([500, 2000])) + s[cut:]
    elif kind < 0.21:  # inversion in the middle
        s = sample(L)
        a1, a2 = len(s) // 3, 2 * len(s) // 3
        s = s[:a1] + revcomp(s[a1:a2]) + s[a2:]
    elif kind < 0.24:  # junk
        s = rand_seq(L)
    else:
        s = sample(L)
    if R.random() < 0.2 and len(s) > 6000:  # diverged (20-30% error) stretch inside the read
        a1 = R.randrange(1000, len(s) - 4000)
        a2 = a1 + R.randrange(1500, 3500)
        s = s[:a1] + mutate(s[a1:a2], R.choice([0.2, 0.3])) + s[a2:]
    if R.random() < 0.2:  # low-identity junk tail (adapter / unalignable end)
        t = rand_seq(R.randrange(300, 3000))
        s = s + t if R.random() < 0.5 else t + s
    if R.random() < 0.5:
        s = revcomp(s)
    s = ont_errors(s, rate=R.choice([0.02, 0.04, 0.07, 0.10, 0.14]))
    reads.append((f"read{i:03d}", s))

os.makedirs(OUT, exist_ok=True)
with gzip.GzipFile(os.path.join(OUT, "ref.fa.gz"), "wb", mtime=0) as f:
    for name, seq in ref.items():
        f.write(f">{name}\n".encode())
        for k in range(0, len(seq), 80):
            f.write((seq[k:k + 80] + "\n").encode())
with gzip.GzipFile(os.path.join(OUT, "reads.fq.gz"), "wb", mtime=0) as f:
    for name, seq in reads:
        f.write(f"@{name}\n{seq}\n+\n{qual(len(seq))}\n".encode())
print(f"{len(ref)} contigs, {sum(map(len, ref.values()))} bp; {len(reads)} reads, "
      f"{sum(len(s) for _, s in reads)} bp -> {OUT}")
