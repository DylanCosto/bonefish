#!/usr/bin/env python3
"""Record-level consistency checks on aligner SAM output (SAM text with header on stdin).

Checks, per record or per read:
  cigar_seq      CIGAR query length (M/I/S/=/X) equals SEQ length
  hardclip       hard clips only on supplementary/secondary records
  one_primary    exactly one primary record per read name
  sa_entries     each SA entry matches another record of the read (contig, position, strand,
                 clips and query/reference spans of the CIGAR, MAPQ, NM), and a read's records
                 list each other
  nm_cs          NM equals the edit distance in the cs tag, and cs agrees with the CIGAR
                 (when --cs was given)
  rg             RG tags refer to an @RG line in the header
  mapq           MAPQ in 0..60 for mapped records; unmapped records carry flag 4 consistently
  mm_ml          MM/ML base-modification tags: on records with the full read sequence (no hard
                 clips), MN equals SEQ length, the MM skip counts fit the read's bases in its
                 original orientation, and ML has one value per MM call; on hard-clipped records
                 MN differs from SEQ length, so tools can tell the tags no longer apply
Prints a table of checks with counts of records examined and failures (first examples shown).
Usage: samtools view -h out.bam | check_sam.py LABEL
"""
import re
import sys
from collections import defaultdict

CIG = re.compile(r"(\d+)([MIDNSHP=X])")
CS = re.compile(r"(:\d+|\*[a-z][a-z]|[+\-][a-z]+|=[A-Za-z]+|~[a-z]{2}\d+[a-z]{2})")
COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def aligned_ops(cigar):
    """(left clip, right clip, query span, reference span): what an SA entry's CIGAR must agree
    on, whether written in full or condensed (minimap2 style, e.g. 100S914M17D20S)."""
    ops = [(int(n), op) for n, op in CIG.findall(cigar)]
    left = ops[0][0] if ops and ops[0][1] in "SH" else 0
    right = ops[-1][0] if len(ops) > 1 and ops[-1][1] in "SH" else 0
    q = sum(n for n, op in ops if op in "MI=X")
    r = sum(n for n, op in ops if op in "MDN=X")
    return (left, right, q, r)


class Checks:
    def __init__(self):
        self.n = defaultdict(int)
        self.bad = defaultdict(int)
        self.ex = defaultdict(list)

    def check(self, name, ok, example=""):
        self.n[name] += 1
        if not ok:
            self.bad[name] += 1
            if len(self.ex[name]) < 3:
                self.ex[name].append(example[:200])


def mm_ok(seq_orig, mm, ml):
    """MM skip counts fit the bases of seq_orig; returns (ok, number of calls)."""
    calls = 0
    for entry in mm.rstrip(";").split(";"):
        if not entry:
            continue
        parts = entry.split(",")
        head = parts[0]
        base = head[0]
        skips = [int(x) for x in parts[1:]]
        calls += len(skips)
        avail = len(seq_orig) if base == "N" else seq_orig.upper().count(base)
        if sum(skips) + len(skips) > avail:
            return False, calls
    return calls == ml, calls


def main():
    label = sys.argv[1] if len(sys.argv) > 1 else "run"
    c = Checks()
    rgs = set()
    reads = defaultdict(list)   # qname -> [(flag, rname, pos, strand, cigar, mapq, nm, sa)]
    for line in sys.stdin:
        if line.startswith("@"):
            if line.startswith("@RG"):
                for t in line.rstrip("\n").split("\t")[1:]:
                    if t.startswith("ID:"):
                        rgs.add(t[3:])
            continue
        f = line.rstrip("\n").split("\t")
        qname, flag, rname, pos, mapq, cigar, seq = f[0], int(f[1]), f[2], int(f[3]), int(f[4]), f[5], f[9]
        tags = {t[:2]: t[5:] for t in f[11:]}
        unmapped = flag & 4
        if unmapped:
            c.check("mapq", cigar == "*" and mapq == 0, line)
            reads[qname].append((flag, None, 0, "", "", 0, 0, ""))
            continue
        c.check("mapq", 0 <= mapq <= 60 and rname != "*" and pos > 0, line)
        ops = CIG.findall(cigar)
        if seq != "*":
            qlen = sum(int(n) for n, op in ops if op in "MIS=X")
            c.check("cigar_seq", qlen == len(seq), f"{qname} {cigar} len={len(seq)}")
        has_h = any(op == "H" for _, op in ops)
        if has_h:
            c.check("hardclip", flag & 0x900 != 0, f"{qname} flag={flag} {cigar}")
        nm = int(tags["NM"]) if "NM" in tags else None
        if "cs" in tags and nm is not None:
            ed = 0
            ref_len = qry_len = 0
            for tok in CS.findall(tags["cs"]):
                if tok[0] == ":":
                    ref_len += int(tok[1:]); qry_len += int(tok[1:])
                elif tok[0] == "=":
                    ref_len += len(tok) - 1; qry_len += len(tok) - 1
                elif tok[0] == "*":
                    ed += 1; ref_len += 1; qry_len += 1
                elif tok[0] == "+":
                    ed += len(tok) - 1; qry_len += len(tok) - 1
                elif tok[0] == "-":
                    ed += len(tok) - 1; ref_len += len(tok) - 1
            cig_ref = sum(int(n) for n, op in ops if op in "MDN=X")
            cig_qry = sum(int(n) for n, op in ops if op in "MI=X")
            c.check("nm_cs", ed == nm and ref_len == cig_ref and qry_len == cig_qry,
                    f"{qname} NM={nm} cs_ed={ed} ref {ref_len}/{cig_ref} qry {qry_len}/{cig_qry}")
        if "RG" in tags:
            c.check("rg", tags["RG"] in rgs, f"{qname} RG={tags['RG']}")
        if "MM" in tags and seq != "*":
            mn = int(tags["MN"]) if "MN" in tags else None
            ml = len(tags["ML"].split(",")) - 1 if "ML" in tags else 0
            if has_h:
                c.check("mm_ml", mn is not None and mn != len(seq), f"{qname} hard-clipped, MN={mn} len={len(seq)}")
            else:
                orig = seq.translate(COMP)[::-1] if flag & 16 else seq
                ok, _ = mm_ok(orig, tags["MM"], ml)
                c.check("mm_ml", ok and (mn is None or mn == len(seq)), f"{qname} MN={mn} len={len(seq)}")
        strand = "-" if flag & 16 else "+"
        reads[qname].append((flag, rname, pos, strand, cigar, mapq, nm, tags.get("SA", "")))

    for qname, recs in reads.items():
        prim = [r for r in recs if not r[0] & 0x900]
        c.check("one_primary", len(prim) == 1, f"{qname}: {len(prim)} primary records")
        linear = [r for r in recs if not r[0] & 0x100 and r[1] is not None]
        keys = {(r[1], r[2], r[3], aligned_ops(r[4]), r[5], r[6]) for r in linear}
        for r in linear:
            sa = [e for e in r[7].split(";") if e]
            if not sa and len(linear) == 1:
                continue
            ok = len(sa) == len(linear) - 1
            for e in sa:
                rn, p, st, cg, mq, enm = e.split(",")
                ok &= (rn, int(p), st, aligned_ops(cg), int(mq), int(enm)) in keys
            c.check("sa_entries", ok, f"{qname} SA={r[7]} records={len(linear)}")

    print(f"## {label}: {len(reads):,} reads")
    print("| check | records examined | failures | examples |")
    print("|---|---:|---:|---|")
    for name in ("cigar_seq", "hardclip", "one_primary", "sa_entries", "nm_cs", "rg", "mapq", "mm_ml"):
        ex = "; ".join(c.ex[name]).replace("|", "/") if c.bad[name] else ""
        print(f"| {name} | {c.n[name]:,} | {c.bad[name]:,} | {ex} |")
    sys.exit(1 if any(c.bad.values()) else 0)


if __name__ == "__main__":
    main()
