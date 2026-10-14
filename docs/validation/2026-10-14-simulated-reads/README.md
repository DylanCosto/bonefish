# Simulated reads with known origin, whole genome (14 Oct 2026)

Real data cannot tell whether an individual read went to the right place. Simulated reads can: each
read's true chromosome, strand and interval is recorded in its name. This checks placement and
mapping quality (MAPQ, which Clair3, Sniffles2 and modkit use to filter reads) for bonefish against
minimap2 2.31 across the whole genome.

## Data and setup

- **Reads:** [Badread](https://github.com/rrwick/Badread) 0.4.2 (nanopore2023 error and quality models)
  from every chromosome of GRCh38 (chr1-22, X, Y, no-alt analysis set), each at the same depth so reads
  fall in proportion to chromosome length (`evidence/sim_wg.sh`). Badread's defaults add about 1% each of
  junk, random and chimeric reads, plus adapters and glitches.

  | set | resembles | identity (mean / max / sd) | length (mean / sd) | depth | reads (Gbp) |
  |---|---|---|---|---:|---:|
  | `r10` | R10.4.1 sup-like | 98 / 99.5 / 1.5% | 15 / 13 kb | 0.33x | 67,308 (1.02) |
  | `lowacc` | older or lower-accuracy reads | 93 / 98 / 4% | 15 / 13 kb | 0.16x | 32,900 (0.49) |
  | `short` | cfDNA / amplicon-like | 98 / 99.5 / 1.5% | 0.5 / 0.3 kb | 0.033x | 196,081 (0.10) |
  | `ul` | ultra-long | 98 / 99.5 / 1.5% | 80 / 50 kb | 0.2x | 7,717 (0.62) |

- **Aligners:** bonefish (this release's code; output identical to v0.1.2) with its prebuilt index;
  minimap2 2.31 `-ax lr:hq`, and also `-ax map-ont` for the low-accuracy set. 16 threads.
- **Scoring** (`evidence/sim_eval.py`): reads with a single true origin are scored (chimeras are
  skipped). A primary alignment is correct if it is on the true chromosome and strand and overlaps the
  true interval. For each MAPQ threshold: reads mapped at or above it and how many of those are wrong
  (the measure of minimap2's `paftools mapeval`). Junk and random reads have no true origin; any of them
  mapped at high MAPQ counts against the aligner.

## Results

| set | aligner | reads scored | unmapped | wrong, all | wrong, MAPQ >= 1 | wrong, MAPQ >= 10 | wrong, MAPQ >= 30 | mapped at MAPQ 60 | junk/random at MAPQ >= 30 |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|
| r10 | **bonefish** | 65,327 | 3,550 | 78 | 20 | **0** | **0** | 60,097 | 35 of 1,331 |
| r10 | minimap2 lr:hq | 65,327 | 3,550 | 78 | 20 | 0 | 0 | 60,096 | 33 of 1,331 |
| ul | **bonefish** | 7,475 | 375 | 2 | 0 | **0** | **0** | 6,967 | 5 of 161 |
| ul | minimap2 lr:hq | 7,475 | 375 | 2 | 0 | 0 | 0 | 6,967 | 5 of 161 |
| short | **bonefish** | 190,294 | 17,864 | 1,606 | 296 | 19 | 1 | 153,224 | 64 of 3,861 |
| short | minimap2 lr:hq | 190,294 | 17,864 | 1,607 | 297 | 19 | 1 | 153,225 | 62 of 3,861 |
| lowacc | **bonefish** | 31,940 | 2,094 | 77 | 44 | 3 | 2 | 28,706 | 33 of 650 |
| lowacc | minimap2 lr:hq | 31,940 | 2,094 | 78 | 45 | 2 | 1 | 28,710 | 33 of 650 |
| lowacc | minimap2 map-ont | 31,940 | 1,837 | 71 | 38 | 0 | 0 | 29,084 | 39 of 650 |

- **R10-like, ultra-long and short reads:** bonefish places reads as minimap2 `lr:hq` does, within one
  or two reads at every MAPQ threshold. On R10-like and ultra-long reads, no read mapped at MAPQ 10 or
  above is in the wrong place, for either aligner.
- **Unmapped reads are mostly reads of nothing:** GRCh38 contains N gaps (centromere and telomere
  gaps, the hard-masked chrY PAR), and Badread simulates reads from them. Reads that are mostly N:
  r10 3,527, lowacc 1,740, short 10,180, ul 393. On the R10-like set, both aligners leave about 23
  real reads unmapped out of about 61,800.
- **Short reads (median ~400 bp):** about 7,700 real reads stay unmapped by both aligners; reads of a
  few hundred bases from repeats give too few unique seeds for a confident chain, and MAPQ is
  correspondingly lower (1,606 wrong at MAPQ 0, 19 at MAPQ >= 10).
- **Low-accuracy reads (~93% identity):** bonefish behaves like minimap2's `lr:hq` preset, which is
  tuned for accurate reads. minimap2's `map-ont` preset (shorter k-mers) maps 257 more of these reads
  and has none wrong at MAPQ >= 5. For low-accuracy data (older chemistries, low-quality basecalls),
  minimap2 `map-ont` remains the better choice; bonefish targets R10.4.1 with current basecall
  models. The real R9.4.1 data in the [variant-calling comparison](../2026-10-14-hac-r9-phasing/README.md)
  (Guppy sup basecalls) shows no loss of accuracy against minimap2 `lr:hq`.
- **Junk and random reads:** both aligners map the same few at high MAPQ (Badread's junk reads are
  low-complexity sequence, which can match genomic tandem repeats).

## Files

- `evidence/results/<set>.<aligner>.md`: the scoring tables; `.curve.tsv`: reads mapped and wrong at
  every MAPQ (for a mapeval-style plot); `.time`: GNU time of each run.
- `evidence/sim_wg.sh` (simulation), `sim_wg_eval.sh` (alignment and scoring), `sim_eval.py`,
  `simulation.log`.
- `evidence/READS_NOT_IN_GIT.sha256`: SHA-256 of the four simulated read sets (1.0 to 3 GB; regenerate
  with `sim_wg.sh`, which fixes Badread's seeds).
- `SHA256SUMS` covers `evidence/`. Local paths are replaced by `${REPO}`, `${DATA}` and `${HOME}`.

Run times in the `.time` files are not a speed benchmark: the runs overlapped with other jobs, and
minimap2 `map-ont` built its index from the FASTA (it was also paused once to free memory).
