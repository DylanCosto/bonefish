#!/usr/bin/env python3
"""Generate the small synthetic spliced dataset used by bonefish/tests/regression.rs (rna_*, cdna_*).

Deterministic (fixed seed). Writes bonefish/tests/data/:
  rna_ref.fa.gz     two contigs with 24 multi-exon genes on both strands (GT..AG introns of
                    80 bp - 20 kb; some genes with alternative isoforms, micro-exons or
                    non-canonical introns), with intergenic sequence and a processed pseudogene
  rna_junc.bed      the annotated transcripts (BED12) for --junc-bed
  rna_reads.fq.gz   direct-RNA-like reads: sense strand, 3' poly(A), 5' truncation, ~6% errors,
                    half of them written with U for T
  cdna_reads.fq.gz  cDNA-like reads: both orientations, ~3% errors, adapter / poly(A) ends

Expected outputs are produced separately (see the test's doc comment).
"""
import gzip
import os
import random

R = random.Random(20261021)
OUT = os.path.join(os.path.dirname(__file__), "..", "bonefish", "tests", "data")
COMP = str.maketrans("ACGT", "TGCA")


def rand_seq(n):
    return "".join(R.choice("ACGT") for _ in range(n))


def revcomp(s):
    return s.translate(COMP)[::-1]


def errors(s, rate):
    """ONT-like errors: substitutions, insertions, deletions (homopolymer-biased)."""
    out = []
    i = 0
    while i < len(s):
        c = s[i]
        x = R.random()
        if x < rate * 0.4:
            out.append(R.choice("ACGT".replace(c, "")))
        elif x < rate * 0.7:
            pass  # deletion
        elif x < rate:
            out.append(c + (c if R.random() < 0.6 else R.choice("ACGT")))
        else:
            out.append(c)
        i += 1
    return "".join(out)


def make_genome():
    contigs = {}
    genes = []  # (contig, strand, [isoform exon lists in genome coordinates], name)
    for ci, name in enumerate(["chrA", "chrB"]):
        seq = [rand_seq(3000)]
        pos = 3000
        for g in range(12):
            gname = f"g{ci}{g:02d}"
            strand = "+" if R.random() < 0.5 else "-"
            n_ex = R.randint(2, 7)
            ex_len = [R.randint(50, 300) for _ in range(n_ex)]
            if n_ex >= 4 and g % 4 == 1:
                ex_len[n_ex // 2] = R.randint(20, 40)  # micro-exon
            in_len = [R.choice([R.randint(80, 400), R.randint(400, 3000), R.randint(3000, 20000)])
                      for _ in range(n_ex - 1)]
            exons = []
            parts = []
            p = pos
            for i in range(n_ex):
                e = rand_seq(ex_len[i])
                exons.append((p, p + ex_len[i]))
                parts.append(e)
                p += ex_len[i]
                if i < n_ex - 1:
                    body = rand_seq(in_len[i] - 4)
                    if g % 6 == 5 and i == 0:
                        donor, acc = "AT", "AC"  # non-canonical (AT-AC) intron
                    else:
                        donor, acc = "GT", "AG"
                    intron = donor + body[:len(body) // 2] + "TTTCTCC" + body[len(body) // 2 + 7:] + acc[:2]
                    intron = intron[:in_len[i] - 2] + acc
                    if strand == "-":
                        # on the genome, a minus-strand intron reads CT..AC (reverse complement of GT..AG)
                        intron = revcomp(intron)
                    parts.append(intron)
                    p += len(intron)
            isoforms = [exons]
            if n_ex >= 4 and g % 3 == 0:
                isoforms.append(exons[:1] + exons[2:])  # exon skipping
            genes.append((name, strand, isoforms, gname))
            seq.append("".join(parts))
            gap = rand_seq(R.randint(2000, 8000))
            seq.append(gap)
            pos = p + len(gap)
        # a processed pseudogene: the spliced copy of this contig's first gene, slightly diverged
        first = [x for x in genes if x[0] == name][0]
        mrna = "".join("".join(seq)[a:b] for a, b in first[2][0])
        seq.append(errors(mrna, 0.03))
        seq.append(rand_seq(3000))
        contigs[name] = "".join(seq)
    return contigs, genes


def transcript(contigs, contig, strand, exons):
    s = "".join(contigs[contig][a:b] for a, b in exons)
    return s if strand == "+" else revcomp(s)


def main():
    contigs, genes = make_genome()
    with gzip.open(os.path.join(OUT, "rna_ref.fa.gz"), "wt", compresslevel=9) as f:
        for name, s in contigs.items():
            f.write(f">{name}\n")
            for i in range(0, len(s), 80):
                f.write(s[i:i + 80] + "\n")
    with open(os.path.join(OUT, "rna_junc.bed"), "w") as f:
        for contig, strand, isoforms, gname in genes:
            for k, ex in enumerate(isoforms):
                start, end = ex[0][0], ex[-1][1]
                sizes = ",".join(str(b - a) for a, b in ex) + ","
                starts = ",".join(str(a - start) for a, b in ex) + ","
                f.write(f"{contig}\t{start}\t{end}\t{gname}.{k}\t0\t{strand}\t{start}\t{end}\t0\t{len(ex)}\t{sizes}\t{starts}\n")
    adapters = ["ATGTACTTCGTTCAGTTACGTATTGCT", "CTTGCCTGTCGCTCTATCTTCAGAGGAG"]
    with gzip.open(os.path.join(OUT, "rna_reads.fq.gz"), "wt", compresslevel=9) as fr, \
         gzip.open(os.path.join(OUT, "cdna_reads.fq.gz"), "wt", compresslevel=9) as fc:
        n = 0
        for contig, strand, isoforms, gname in genes:
            for ex in isoforms:
                tx = transcript(contigs, contig, strand, ex)
                for rep in range(3):
                    n += 1
                    # direct RNA: sense strand, 5' truncated, poly(A) tail
                    cut = R.randint(0, len(tx) // 3)
                    rd = errors(tx[cut:] + "A" * R.randint(10, 40), 0.06)
                    if n % 2 == 0:
                        rd = rd.replace("T", "U")
                    fr.write(f"@rna{n:03d}_{gname}\n{rd}\n+\n{'5' * len(rd)}\n")
                    # cDNA: full length with adapters, either orientation
                    cd = adapters[0] + errors(tx, 0.03) + "A" * R.randint(15, 30) + adapters[1]
                    if R.random() < 0.5:
                        cd = revcomp(cd)
                    fc.write(f"@cdna{n:03d}_{gname}\n{cd}\n+\n{'9' * len(cd)}\n")
        # unmappable junk
        for j in range(4):
            s = rand_seq(R.randint(300, 1500))
            fr.write(f"@rna_junk{j}\n{s}\n+\n{'5' * len(s)}\n")
            fc.write(f"@cdna_junk{j}\n{s}\n+\n{'9' * len(s)}\n")


if __name__ == "__main__":
    main()
