# Other references: T2T-CHM13 and GRCh38 with ALT contigs (14 Oct 2026)

bonefish against minimap2 2.31 (`-ax lr:hq`) on two references besides the GRCh38 no-alt analysis set
used elsewhere, with the genome-wide HG001 sample (161,762 R10.4.1 sup reads, 2.71 Gbp), 16 threads,
the machine of the [aligner comparison](../2026-10-14-aligner-comparison/README.md).

- **T2T-CHM13v2.0** (UCSC `hs1.fa.gz`, 25 sequences, 3.12 Gbp).
- **GRCh38 with ALT contigs and decoys** (NCBI `GCA_000001405.15_GRCh38_full_plus_hs38d1_analysis_set`,
  2,841 sequences: the primary assembly, 261 ALT contigs, 2,385 decoys, EBV). The 261 contig names ending
  in `_alt` are given to both aligners with `--alt` (`evidence/GRCh38_alt.list`), which lowers the
  mapping quality of hits that compete only with ALT copies.

## Results

| reference | aligner | wall | CPU time | peak RSS | reads mapped | bases aligned | identity | MAPQ >= 20 |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| T2T-CHM13 | bonefish | 41 s | 587 s | 14.8 GB | 96.01% | 97.33% | 0.9823 | 94.42% |
| | minimap2 | 166 s | 2,136 s | 18.6 GB | 96.01% | 97.46% | 0.9817 | 94.40% |
| GRCh38 + ALT + decoys (`--alt`) | bonefish | 49 s | 704 s | 15.2 GB | 96.35% | 96.52% | 0.9790 | 95.11% |
| | minimap2 | 206 s | 2,568 s | 19.2 GB | 96.35% | 96.80% | 0.9773 | 94.96% |
| GRCh38 no-alt (for reference) | bonefish | 46 s | 658 s | 15.0 GB | 96.30% | 96.44% | 0.9791 | 94.96% |
| | minimap2 | 192 s | 2,563 s | 18.3 GB | 96.30% | 96.78% | 0.9773 | 94.84% |

ALT handling (`evidence/alt_stats.tsv`):

| | primaries on ALT contigs | on decoys | primaries in the MHC (chr6:28.5-33.4 Mb) | of which MAPQ >= 20 |
|---|---:|---:|---:|---:|
| bonefish | 155 | 1,077 | 217 | 213 (98.2%) |
| minimap2 | 154 | 1,073 | 217 | 213 (98.2%) |

Index build (wall): CHM13 bonefish 49 s, minimap2 57 s; GRCh38 + ALT bonefish 84 s, minimap2 58 s.

- bonefish is 4.0x (CHM13) and 4.2x (with ALT contigs) faster than minimap2, as on the no-alt reference.
- Both map the same share of reads on each reference. On CHM13 alignment identity is higher for both
  (the reads match the more complete reference better).
- With ALT contigs and decoys present, bonefish and minimap2 place the same number of primaries on ALT
  contigs and decoys, and keep the same mapping quality in the MHC, where most ALT contigs lie.

The mapping columns are not variant-calling accuracy; that was measured on GRCh38 no-alt
([GIAB whole genomes](../2026-10-12-giab-wgs/README.md)).

## Files

`evidence/stats/<run>/bam_summary.tsv` (mapping statistics), `evidence/runs/*.time` (GNU time of index
builds and mapping), `*.err`, `run.log`, `alt_stats.tsv`, `GRCh38_alt.list`, `ref_tests.sh`, and
`REFERENCES.sha256` (the downloaded reference files). `SHA256SUMS` covers `evidence/`. Local paths are
replaced by `${REPO}`, `${DATA}` and `${HOME}`.
