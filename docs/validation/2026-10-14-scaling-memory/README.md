# Thread scaling and memory (14 Oct 2026)

How bonefish's speed changes with thread count, and how little memory a human genome can be
aligned in, for laptops and shared servers.

## Setup

- Reads: the genome-wide HG001 sample (161,762 R10.4.1 sup reads, 2.7 Gbp); reference GRCh38 no-alt.
- Machine: Intel Core i5-12600K (6 performance cores with 2 threads each, 4 efficiency cores: 16
  threads), 31 GB RAM, Linux. Nothing else running (the script waits for an idle machine).
- SAM output is hashed (MD5 without `@PG`) instead of written, so disk speed does not count.
  Index files were read into the page cache first. Script: `evidence/scaling.sh`.
- bonefish with its raw index (memory-mapped), minimap2 2.31 `-ax lr:hq` with its `.mmi` index.

## Threads

| threads | bonefish wall | bonefish CPU | bonefish peak memory | minimap2 wall | minimap2 CPU | minimap2 peak memory |
|---:|---:|---:|---:|---:|---:|---:|
| 1 | 360 s | 372 s | 13.5 GB | | | |
| 2 | 193 s | 398 s | 13.6 GB | | | |
| 4 | 103 s | 424 s | 13.8 GB | 385 s | 1,503 s | 14.3 GB |
| 8 | 68 s | 542 s | 14.3 GB | 234 s | 1,746 s | 15.8 GB |
| 12 | 55 s | 627 s | 14.7 GB | | | |
| 16 | 51 s | 667 s | 15.0 GB | 191 s | 2,569 s | 18.7 GB |

- The output is identical at every thread count (the same MD5 as the published v0.1.1 output,
  `deed7949d35f…`).
- bonefish is 3.4-3.7x faster than minimap2 at the same thread count, and on 4 threads it is faster
  than minimap2 on 16.
- Speed-up from 1 to 16 threads is 7x. CPU time per thread rises past 4 threads, as expected on this
  CPU: threads beyond the 6 performance cores share cores (hyper-threading) or run on the slower
  efficiency cores.

## Memory

Peak memory has two parts: the index (the raw GRCh38 index is 8.5 GB, memory-mapped) and the batch
of reads being aligned. The batch size is `-K` (query bases per batch).

| setting (4 threads unless noted) | wall | peak memory | output |
|---|---:|---:|---|
| raw index, `-K` 500M (the old default) | 103 s | 13.8 GB | identical |
| `.idx` index, `-K` 500M | 110 s | 15.1 GB | identical |
| raw index, `-K 10M` | 107 s | **8.6 GB** | identical |
| `.idx` index, `-K 10M` | 107 s | 9.1 GB | identical |
| index split into 2 parts (`-I 2G`, built from the FASTA each run) | 893 s | 7.5 GB | differs |
| index split into 4 parts (`-I 1G`) | 2,114 s | 4.3 GB | differs |

At 16 threads (`evidence/batch_size.tsv`):

| `-K` | SAM wall | SAM peak | sorted BAM wall | sorted BAM peak |
|---|---:|---:|---:|---:|
| 500M (old default) | 49.4 s | 15.0 GB | 88.0 s | 16.3 GB |
| **100M (new default)** | 48.3 s | 12.3 GB | 84.1 s | 14.9 GB |
| 50M | 50.1 s | 12.1 GB | | |
| 20M | 51.0 s | 10.8 GB | | |
| 10M | 52.5 s | 10.4 GB | | |

- **The default batch is now 100 Mbp** (was 500 Mbp): the same speed with 2.7 GB less memory for SAM
  output, 1.5 GB less for sorted BAM. The batch size does not change the output.
- **A human genome on a 16 GB machine:** `bonefish -t 4 -K 10M` peaks at 8.6 GB and maps this 2.7 Gbp
  sample in 107 s, so a 30x genome (about 93 Gbp) would take about an hour on 4 threads, with
  the same output as on a 16-thread server. (SAM output; a sorted BAM adds the sort buffers, which
  `--sort-mem` limits.)
- **Split index (`-I`):** fits in 4-8 GB but is 9-20x slower (each part is built from the FASTA and
  every read is aligned against every part) and, as with minimap2, changes the output slightly
  (mapping quality is computed per part). Use it only when memory is under 8 GB.

## Files

`evidence/scaling.tsv` (all runs), `batch_size.tsv`, `scaling.sh`, and the last lines of each run's log
in `evidence/logs/`. `SHA256SUMS` covers `evidence/`.
