#!/usr/bin/env python3
"""Extract the small HiFi dataset used by bonefish/tests/regression.rs (hifi_*).

Eight HG002 PacBio Revio HiFi reads (PacBio's public 2022Q4 HG002 rep1 run,
m84011_220902_175841_s1.hifi_reads.bam) across three GIAB v5.0q structural variants, where
lr:hq scoring splits the variant into several false deletions and -x hifi does not:
  chr1:187495694  DEL 1900 + INS 248 (the dna preset: six deletions of 87-441 bp)
  chr4:61133657   DEL 2575           (two deletions, 1333 + 1226)
  chr4:129144970  INS 317 + DEL 2315 (four deletions and ~300 bp of junk alignment)
with 27-30 kb of GRCh38 around each.

Usage: make_regression_hifi_data.py READS.fq[.gz] GRCh38.fa
  READS: FASTQ holding the reads below in sequencing orientation (e.g. `samtools fastq` of the
  run's BAM). GRCh38.fa: GRCh38 (no-alt analysis set), indexed with `samtools faidx`.
Writes bonefish/tests/data/hifi_ref.fa.gz (one contig per region, named chrN_START) and
hifi_reads.fa.gz (FASTA: HiFi base qualities compress poorly and are not used by the aligner).
"""
import gzip
import os
import subprocess
import sys

RUN = "m84011_220902_175841_s1"
READS = ["176294863", "25235440",                 # chr1:187495694
         "155260449", "78385068", "220463385",    # chr4:61133657
         "50401584", "175573594", "120133061"]    # chr4:129144970
REGIONS = [("chr1", 187485000, 187512000), ("chr4", 61120000, 61147000), ("chr4", 129130000, 129160000)]
OUT = os.path.join(os.path.dirname(__file__), "..", "bonefish", "tests", "data")


def wrap(s, n=80):
    return "\n".join(s[i:i + n] for i in range(0, len(s), n))


def main():
    reads_path, ref = sys.argv[1], sys.argv[2]
    names = {f"{RUN}/{r}/ccs" for r in READS}
    found = {}
    opener = gzip.open if reads_path.endswith(".gz") else open
    with opener(reads_path, "rt") as f:
        while len(found) < len(names):
            h = f.readline()
            if not h:
                break
            s, _, _ = f.readline().strip(), f.readline(), f.readline()
            name = h[1:].split()[0]
            if name in names:
                found[name] = s
    missing = names - found.keys()
    if missing:
        sys.exit(f"reads not found: {sorted(missing)}")
    # gzip mtime 0: byte-identical output on every run
    with gzip.GzipFile(os.path.join(OUT, "hifi_reads.fa.gz"), "wb", mtime=0) as g:
        for r in READS:
            n = f"{RUN}/{r}/ccs"
            g.write(f">{n}\n{wrap(found[n])}\n".encode())
    with gzip.GzipFile(os.path.join(OUT, "hifi_ref.fa.gz"), "wb", mtime=0) as g:
        for c, s, e in REGIONS:
            fa = subprocess.run(["samtools", "faidx", ref, f"{c}:{s + 1}-{e}"],
                                capture_output=True, text=True, check=True).stdout
            seq = "".join(fa.splitlines()[1:])
            g.write(f">{c}_{s}\n{wrap(seq)}\n".encode())


if __name__ == "__main__":
    main()
