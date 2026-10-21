#!/usr/bin/env python3
"""Score transcriptome alignments of Badread reads simulated from GENCODE transcripts.

Read names carry the source transcript (`ENST...|type|GENE,+strand,a-b`); reference names are
GENCODE transcriptome headers (`ENST...|ENSG...|...|...|tx_name|GENE|len|type|`). Per BAM: reads
mapped, primary on the true transcript, primary on a transcript of the true gene, and the true
transcript among any of the read's alignments (what EM quantifiers use).

Usage: tx_sim_truth.py READS.fq.gz A.bam [B.bam ...]
"""
import gzip
import subprocess
import sys


def main():
    truth = {}
    with gzip.open(sys.argv[1], "rt") as fq:
        for i, line in enumerate(fq):
            if i % 4 == 0:
                rid, _, desc = line[1:].rstrip("\n").partition(" ")
                src = desc.split(",")[0].split("|")
                truth[rid] = (src[0], src[2] if len(src) > 2 else "")
    n = len(truth)
    for bam in sys.argv[2:]:
        prim = {}
        anyhit = set()
        p = subprocess.Popen(["samtools", "view", bam], stdout=subprocess.PIPE, text=True)
        for line in p.stdout:
            f = line.split("\t", 4)
            flag = int(f[1])
            if flag & 4 or f[0] not in truth:
                continue
            ref = f[2].split("|")
            tx, gene = ref[0], ref[5] if len(ref) > 5 else ""
            if tx == truth[f[0]][0]:
                anyhit.add(f[0])
            if not flag & 0x900:
                prim[f[0]] = (tx, gene)
        p.wait()
        mapped = len(prim)
        tx_ok = sum(prim[r][0] == truth[r][0] for r in prim)
        gene_ok = sum(prim[r][1] == truth[r][1] for r in prim)
        print(f"{bam}: mapped {100 * mapped / n:.2f}%  primary on true transcript {100 * tx_ok / n:.2f}%  "
              f"on true gene {100 * gene_ok / n:.2f}%  true transcript among alignments {100 * len(anyhit) / n:.2f}%")


if __name__ == "__main__":
    main()
