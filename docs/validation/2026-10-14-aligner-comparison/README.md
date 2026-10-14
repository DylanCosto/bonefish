# bonefish vs other ONT aligners (14 Oct 2026)

Speed, memory and mapping statistics of eight aligners on the same R10.4.1 reads, reference and
machine. The tables are summarised in [docs/benchmarks.md](../../benchmarks.md); this folder holds
the raw measurements.

## Setup

- **Machine:** Intel Core i5-12600K (16 threads, AVX2, no AVX-512), 31 GB RAM, Samsung 980 NVMe,
  Linux Mint (kernel 7.0). Desktop applications closed and swap reset before the runs.
- **Reference:** GRCh38 no-alt analysis set (195 sequences, 3.10 Gbp).
- **Reads:** HG001, R10.4.1 dorado sup, a random genome-wide sample of 161,762 reads (2.71 Gbp,
  N50 28 kb; `samtools view --subsample 0.027 --subsample-seed 7` of ONT's GIAB BAM, as FASTQ).
  The 10% subset is every tenth read of it (16,177 reads, 270 Mbp).
- **Threads:** 16 for every tool. One run per tool and set. GNU `time -v` around the aligner.
- **Output:** SAM to a pipe that keeps columns 1-6 and the tags (`cut -f1-6,12-`), so disk writes
  do not dominate; statistics from those columns with `bam_summary.py`. The sorted-BAM runs write
  BAM + index to the same NVMe disk (NTFS through ntfs-3g, ~300 MB/s).

| tool | version | preset / options | index |
|---|---|---|---|
| bonefish | 0.1.1 (+ input fixes, identical output) | `dna` (default) | raw index (`-d X.raw --idx-raw`) |
| minimap2 | 2.31-r1302 | `-ax lr:hq` | `.mmi` |
| dorado aligner | 2.1.2 | default (`lr:hq`), `--emit-sam`; `-o DIR` for sorted BAM | minimap2's `.mmi` |
| mm2plus | 1.3 (bioconda) | `-ax lr:hq` | own `.mmi` |
| mm2-fast | git `bwa-mem2/mm2-fast` (minimap2 2.24), AVX2 build | `-ax map-ont` (no `lr:hq` in 2.24) | own `.mmi` |
| Winnowmap2 | 2.03 (bioconda) | `-W rep15.txt -ax map-ont` | built per run (see notes) |
| NGMLR | 0.2.7 (bioconda) | `-x ont` | `.ngm` files |
| LRA | 1.3.7.1 (bioconda) | `align -ONT -p s` | `lra index -ONT` |

Exact command lines are in `evidence/aligner_bench.sh`; versions in `evidence/versions.txt`.

## Results

**Genome-wide sample, SAM output**

| aligner | wall | CPU time | peak RSS | reads mapped | bases aligned | identity | MAPQ >= 20 | supplementary |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| bonefish | 46.0 s | 658 s | 15.0 GB | 96.30% | 96.44% | 0.9791 | 94.96% | 14,632 |
| minimap2 | 191.7 s | 2,563 s | 18.3 GB | 96.30% | 96.78% | 0.9773 | 94.84% | 13,543 |
| dorado aligner | 197.7 s | 2,167 s | 13.3 GB | 96.30% | 95.85% | 0.9791 | 95.09% | 22,430 |
| mm2plus | 236.3 s | 3,059 s | 18.5 GB | 96.30% | 96.78% | 0.9773 | 94.84% | 13,554 |
| mm2-fast | 641.3 s | 9,488 s | 19.6 GB | 96.99% | 97.14% | 0.9741 | 94.50% | 11,891 |
| Winnowmap2 | 1,764.4 s | 25,356 s | 19.4 GB | 96.19% | 96.83% | 0.9797 | 95.12% | 21,334 |

**10% subset, SAM output**

| aligner | wall | CPU time | peak RSS | reads mapped | bases aligned | identity |
|---|---:|---:|---:|---:|---:|---:|
| bonefish | 11.1 s | 68 s | 10.1 GB | 96.30% | 96.57% | 0.9788 |
| dorado aligner | 27.4 s | 210 s | 12.8 GB | 96.30% | 96.01% | 0.9789 |
| mm2plus | 33.6 s | 299 s | 16.8 GB | 96.30% | 96.92% | 0.9770 |
| minimap2 | 34.3 s | 268 s | 16.1 GB | 96.30% | 96.92% | 0.9770 |
| mm2-fast | 67.7 s | 932 s | 17.1 GB | 96.92% | 97.21% | 0.9740 |
| LRA | 74.6 s | 702 s | 11.3 GB | 89.56% | 94.71% | 0.9896 |
| Winnowmap2 | 202.8 s | 2,307 s | 14.2 GB | 96.08% | 96.94% | 0.9794 |
| NGMLR | 2,908.9 s | 24,561 s | 12.8 GB | 93.00% | 95.59% | 0.9855 |

**Reads to sorted, indexed BAM (genome-wide sample)**

| pipeline | wall | CPU time | peak RSS | records |
|---|---:|---:|---:|---:|
| `bonefish -t 16 -o out.bam` | 54.8 s | 758 s | 16.4 GB | 212,292 |
| `minimap2 -ax lr:hq \| samtools sort -@16 -m 384M --write-index` | 218.5 s | 2,699 s | 18.6 GB | 214,121 |
| `dorado aligner -t 16 -o DIR` | 573.5 s | 2,501 s | 13.3 GB | 222,574 |

All three BAMs hold the 161,762 input reads as primary records, are coordinate-sorted and indexed.

**Index build** (wall / peak RSS): bonefish 61.8 s / 10.5 GB; minimap2 64.7 s / 11.3 GB; mm2plus
61.4 s / 13.3 GB; mm2-fast 63.4 s / 11.3 GB; Winnowmap2 63.3 s / 4.6 GB plus meryl k-mer counting
130.3 s / 3.2 GB; LRA 196.0 s / 16.2 GB; NGMLR 416.4 s / 9.5 GB.

## Notes and limits

- These are speed and mapping statistics, not accuracy. bonefish's accuracy has been measured
  with variant calling against GIAB truth sets, against minimap2 only
  ([GIAB whole genomes](../2026-10-12-giab-wgs/README.md)).
- "Identity" is 1 - NM / alignment columns over primary alignments: tools that clip more of each
  read score higher. "Bases aligned" is the share of primary-read bases inside the alignment.
- Winnowmap 2.03 produced no alignments from a prebuilt index on test data, so it indexed the
  FASTA on every run; `idx_winnowmap` (a run on 10 reads) measures that cost (63 s).
- mm2-fast's speedups target AVX-512; on this AVX2 CPU, with the older `map-ont` preset (wider
  extension bands than `lr:hq`), it was slower than minimap2.
- dorado aligner reports far more supplementary alignments than minimap2 on the same index, and
  aligns fewer read bases, so its default options differ from minimap2's `-x lr:hq`. Its
  sorted-output mode ran at ~32% CPU on average.
- Times on the 10% subset include loading the index, which favours tools with small or quickly
  loaded indexes on short runs; compare CPU time there.
- Swap rose during the genome-wide minimap2 run (2.5 to 7.9 GB) and the dorado sorted-BAM run (to
  8.4 GB); see `evidence/resources.csv` (10 s samples, `stage` column). minimap2's time agrees
  with earlier measurements on this machine (201.5 s).

## Files

- `evidence/versions.txt`, `evidence/driver.log`: tool versions and the run log.
- `evidence/runs/<run>.time`: GNU `time -v` of each run; `<run>.err`: the tool's stderr.
- `evidence/stats/<run>/bam_summary.tsv`: mapping statistics.
- `evidence/resources.csv`: CPU, memory, swap and I/O every 10 s.
- `evidence/aligner_bench.sh`, `monitor.sh`, `bam_summary.py`: the scripts (paths are those of the
  development machine).
- `evidence/OUTPUTS_NOT_IN_GIT.sha256`: SHA-256 of the BAMs and per-read tables kept outside Git.
- `SHA256SUMS`: checksums of everything in `evidence/` (`sha256sum -c SHA256SUMS`).

Local paths in the logs are replaced by `${REPO}` (the development repository), `${DATA}` (the
data volume), `${HOME}` and `${SCRATCH}`.
