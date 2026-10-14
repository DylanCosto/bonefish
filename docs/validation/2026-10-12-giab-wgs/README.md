# Whole-genome validation on GIAB HG002 and HG005 (12 Oct 2026)

bonefish against minimap2 2.31, from reads to a sorted, indexed BAM and on to small-variant and
structural-variant calls scored against GIAB truth sets. Summarised in
[docs/benchmarks.md](../../benchmarks.md); this folder holds the per-sample reports and raw
results. The bonefish build used gives output identical to v0.1.1 on x86 (and, since v0.1.1, on
arm64).

## Data and methods

- **Reads:** ONT's GIAB 2025.01 release, R10.4.1, dorado sup v5.0.0. HG005: all reads of flow cell
  PAW87816 (27.9x). HG002: a 60% subsample of PAW71238 (27.2x; `samtools view -s 7.6`). Primary and
  unmapped records of ONT's BAM as FASTQ, the same file for both aligners.
- **Alignment (16 threads, `-K 50M`):** `minimap2 -ax lr:hq | samtools sort -@16 -m 384M --write-index`
  against `bonefish -a --sort` (now simply `bonefish -o out.bam`).
- **Small variants:** Clair3 2.0.3, model `r1041_e82_400bps_sup_v500`, longphase phasing;
  `rtg vcfeval` against GIAB v4.2.1 (overall and GIAB v3.3 stratifications), and for HG002 also
  against GIAB v5.0q (the T2T-Q100-based benchmark, 2.74 Gbp including chrX/Y).
- **Structural variants (HG002):** Sniffles2, Truvari against GIAB v5.0q SVs, strict (one-to-one)
  and multi-match.
- **Machine:** as in the [aligner comparison](../2026-10-14-aligner-comparison/README.md); CPU,
  memory, swap and I/O sampled every 10 s (`resources.csv`).

## Results

| | HG005 minimap2 | HG005 bonefish | HG002 minimap2 | HG002 bonefish |
|---|---:|---:|---:|---:|
| reads -> sorted, indexed BAM | 114.7 min | 34.8 min | 114.8 min | 37.5 min |
| CPU time (aligner + sort) | 25.3 h | 7.2 h | 25.3 h | 7.2 h |
| peak RSS | 16.7 GB + sort 7.3 GB | 19.6 GB | 16.9 GB + sort 7.2 GB | 19.7 GB |
| SNP F1 (TP / FP / FN) | 0.9975 (3,268,693 / 8,774 / 7,778) | 0.9974 (3,268,656 / 8,915 / 7,815) | 0.9971 (3,358,537 / 11,752 / 7,837) | 0.9971 (3,358,505 / 11,851 / 7,869) |
| indel F1 | 0.9395 | 0.9396 | 0.8737 | 0.8737 |
| segmental-duplication SNP F1 (FP) | 0.9735 (4,594) | 0.9728 (4,750) | 0.9618 (7,067) | 0.9610 (7,227) |
| SV F1 vs v5.0q, strict / multi-match | – | – | 0.8055 / 0.8706 | 0.8094 / 0.8741 |
| primary alignments identical / moved | | 90.6% / 1.1% | | 94.9% / 1.4% |
| reads newly unmapped | | 3 | | 2 |
| Clair3 wall (CPU) | 74.7 min (14.4 h) | 81.3 min (15.4 h) | 80.4 min (15.4 h) | 89.7 min (17.0 h) |

HG002 small variants against GIAB v5.0q:

| region | SNP F1 minimap2 (FP) | SNP F1 bonefish (FP) | indel F1 minimap2 | indel F1 bonefish |
|---|---:|---:|---:|---:|
| all | 0.9862 (46,658) | 0.9863 (45,136) | 0.7075 | 0.7076 |
| new vs v4.2.1 | 0.8774 (42,180) | 0.8787 (40,782) | 0.4475 | 0.4475 |
| segmental duplications | 0.9333 (8,857) | 0.9337 (8,751) | 0.6506 | 0.6512 |
| low mappability | 0.9660 | 0.9662 | 0.6488 | 0.6493 |
| tandem repeats / homopolymers | 0.8942 | 0.8947 | 0.6222 | 0.6222 |
| not in difficult regions | 0.9971 | 0.9972 | 0.9584 | 0.9586 |

SNP and indel F1 are equal genome-wide on both samples (|ΔF1| <= 0.0001). Against v4.2.1
bonefish has about 160 more SNP false positives in segmental duplications per genome; against
v5.0q it has fewer false positives there and overall. SV recall is higher on HG002 (+262 TP,
+140 FP). Clair3 takes 9-12% longer on bonefish BAMs. bonefish aligns slightly fewer read bases
(96.2% vs 96.6%) at higher identity (0.979 vs 0.977): low-identity read ends stay soft-clipped.

## Files

Per sample (`evidence/HG002`, `evidence/HG005`): `report.md` (the generated report), `run.log`,
`resources.csv`, `compare_primary.tsv` (primary alignments, bonefish vs minimap2),
`sv_agreement.tsv`, and per aligner: `accuracy.tsv`, `strat.tsv`, `bam_summary.tsv`,
`clair3_counts.tsv`, `flagstat.txt`, `sv_accuracy.tsv`, rtg (`vcfeval_*.summary.txt`) and Truvari
(`truvari_*.summary.json`) summaries, GNU time of every stage; HG002 also `q100/` (v5.0q rtg
summaries). `OUTPUTS_NOT_IN_GIT.sha256` lists the SHA-256 of the VCFs and per-read tables kept
outside Git. Scripts: `sample_eval.sh`, `sample_report.py`, `q100_eval.sh`, `compare_tables.py`.
`SHA256SUMS` covers `evidence/`. Local paths are replaced by `${REPO}`, `${DATA}`, `${HOME}`.
