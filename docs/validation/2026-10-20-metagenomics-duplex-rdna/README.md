# Metagenomics, duplex reads and rDNA (20 Oct 2026)

Three kinds of ONT DNA data not covered by the earlier human whole-genome validations, each mapped
with bonefish v0.1.7 (release binary) and minimap2 2.31 (`-ax lr:hq`), 16 threads, the same reads.
There is no truth set for these reads, so the comparison is against minimap2: placement, mapping
quality, identity, abundance estimates and speed.

## Data

- **Metagenome:** ONT's Zymo fecal reference (`s3://ont-open-data/zymo_fecal_2025.05`), flow cell
  PAU85136, dorado sup 5.0.0 with 4mC/5mC/6mA calls; the first 12 GB of the unaligned BAM, 1,701,552
  reads (10.0 Gbp, mean 5.9 kb), converted to FASTQ with the `MM`/`ML` tags. Reference: ONT's
  metaMDBG assembly of eight flow cells of the same sample (`analysis/EigthFC/bgzip_contigs.fa.gz`):
  29,984 contigs, 1.10 Gbp, 226 circular; contig-to-bin table from DAS Tool (380 bins with reads).
- **Duplex:** HG002 R10.4.1 stereo duplex reads from UCSC (HPRC submission
  `UCSC_HG002_R1041_Duplex_Dorado`, Dorado 0.1.1, flow cell 3 of 11/15/22, Q10 pass): 356,457 reads,
  9.4 Gbp, mean 26 kb. Reference: GRCh38 no-alt analysis set.
- **rDNA:** reads of the genome-wide HG001 sample (`hg001_wg_1x.fq.gz`) with any alignment to
  `chrUn_GL000220v1`, which holds GRCh38's 45S rDNA unit, by either aligner: 530 reads. Mapped to
  GRCh38 and to T2T-CHM13 (hs1), where the rDNA arrays are modelled on the short arms of chr13, 14,
  15, 21 and 22.

## Results

| | metagenome: minimap2 | bonefish | duplex: minimap2 | bonefish |
|---|---:|---:|---:|---:|
| reads mapped | 97.86% | 97.86% | 98.43% | 98.43% |
| MAPQ >= 20 | 86.78% | 86.79% | 93.43% | 93.66% |
| primaries identical to minimap2 | | 99.75% | | 96.33% |
| primaries at an overlapping locus | | 99.95% | | 98.37% |
| identity (BLAST-like / gap-compressed) | 0.9832 / 0.9921 | 0.9834 / 0.9922 | 0.9855 / 0.9939 | 0.9874 / 0.9942 |
| aligned bases, primary + supplementary | 9.063 Gb | 9.059 Gb | 9.498 Gb | 9.453 Gb |
| wall time, SAM output (`walltime.txt`) | 372 s | **169 s** | 715 s* | **186 s** |
| CPU time | 5,133 s | **2,008 s** | 9,657 s | **2,369 s** |
| peak memory | 13.0 GB | 9.7 GB | 19.9 GB | 14.7 GB |

\* minimap2's duplex run wrote BAM through `samtools view -b`; at 13.5 cores on average it was not
limited by that step. The metagenome runs through `samtools view -b` were (745 s and 769 s for 17 GB
of BAM with modified-base tags), so the wall times above are from runs writing SAM to `/dev/null`.

- **Metagenome:** 99.75% of primaries are identical to minimap2's. Per-contig read counts correlate
  with Pearson r = 1.000000 (0.99972 on log counts); 555 reads (0.033%) would change contig and 199
  (0.012%) bin, and no bin with 100 or more reads changes by more than 3%
  (`metagenomics/meta.abundance.txt`). Reads spanning the origin of circular contigs are split into
  primary and supplementary alignments as with minimap2 (8,032 against 8,015). Every mapped read
  keeps its `MM`/`ML` tags with `-y`. bonefish needs 2.6x less CPU time; the gain is smaller than on
  human genomes (about 4x) because these reads are shorter and less of the time is in long gaps.
- **Duplex:** mapped share equal, mapping quality slightly higher, primary identity higher (low-identity
  read ends are left soft-clipped, as for simplex reads: 0.5% fewer aligned bases). Of the 1.6% of
  primaries placed differently, 80% have MAPQ below 5 in both aligners (repeat copies); 332 reads
  (0.09%) are placed confidently but differently by the two (`duplex/moved_reads.txt`). 4.1x less CPU
  time, 3.8x faster.
- **rDNA:** all 530 reads are mapped by both aligners on GRCh38; on CHM13 both leave the same 39
  short reads (280-560 bp) unmapped. 98.5% (GRCh38) and 98.9% (CHM13) of primaries are at the same
  place in both, with the same mapping-quality profile: on CHM13, 224 reads have MAPQ 0 in both,
  as expected for hundreds of near-identical 43 kb copies (`rdna/`). For rDNA studies, aligning to
  a single rDNA unit remains the usual approach, for either aligner.

## Files

`evidence/metagenomics/` (statistics, primary comparison, abundance comparison, timings, circular
contig list, contig-to-bin table), `evidence/duplex/` (statistics, comparison, moved-read MAPQ,
timings), `evidence/rdna/` (summary and per-read agreement), `evidence/walltime.txt`,
`evidence/versions.txt`, and `evidence/scripts/` (`meta_duplex_fetch.sh` downloads the data,
`meta_duplex_eval.sh` runs both aligners and the statistics, `walltime.sh`, `rdna_check.sh`,
`aln_stats.py`, `abundance_compare.py`, `compare_primary.py`, `s3_parallel_cat.py`). Reference
sha256: Zymo fecal contigs `f5b108a3ca127862828258a8126b73e12dfe876da1f6b7d798de4609120c5c39`.
`SHA256SUMS` covers `evidence/`. Local paths are replaced by `${BENCH}`, `${DATA}`, `${REPO}` and
`${HOME}`.
