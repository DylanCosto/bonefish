# CpG methylation from bonefish vs minimap2 BAMs (14 Oct 2026)

Does the aligner change the methylation calls? HG002 reads basecalled with modified-base models
were aligned by minimap2 and by bonefish, and modkit tabulated CpG 5mC/5hmC from each BAM.

## Data and methods

- **Reads:** every primary read on chr20 of ONT's HG002 flow cell PAW70337 (GIAB 2025.01, R10.4.1,
  dorado sup v5.0.0 with 5mCG/5hmCG and 6mA models): 188,140 reads, 3.52 Gbp (~55x of chr20),
  each with `MM`/`ML`/`MN` tags. The same FASTQ (`samtools fastq -T MM,ML,MN`) went to both aligners,
  which mapped it to the whole GRCh38 no-alt genome.
- **minimap2 2.31:** `minimap2 -t 16 -ax lr:hq -y | samtools sort --write-index`.
- **bonefish:** `bonefish -t 16 -y -o out.bam` (sorted BAM and index in one step); also run on the
  original BAM directly (`bonefish -y -o out.bam ref reads.bam`) to check the BAM-input path.
- **modkit 0.6.4:** `modkit pileup --modified-bases 5mC 5hmC --cpg --combine-strands --region chr20
  --ref GRCh38`, default filtering. modkit's automatic call threshold came out the same for all three
  BAMs (0.7285 for C), so the calls are filtered identically.
- **Comparison** (`modkit_compare.py`): per CpG site with valid coverage >= N in both pileups.

## Results

| valid coverage >= | CpG sites (minimap2 / bonefish / both) | 5mC Pearson r | mean abs. diff (points) | sites differing > 10 / > 25 points | 5hmC Pearson r | 1 kb bins r (5mC) |
|---:|---|---:|---:|---:|---:|---:|
| 5 | 760,925 / 760,652 / 759,089 | 0.9981 | 0.29 | 0.98% / 0.18% | 0.9845 | 0.9991 |
| 10 | 752,852 / 753,657 / 751,483 | 0.9988 | 0.21 | 0.66% / 0.10% | 0.9924 | 0.9991 |
| 20 | 731,467 / 731,850 / 730,056 | 0.9994 | 0.13 | 0.36% / 0.04% | 0.9970 | 0.9991 |

- chr20 CpG 5mC level: 64.38% from minimap2, 64.42% from bonefish; 5hmC 1.06% from both.
  Mean valid coverage per CpG 45.3 vs 45.1 (r = 0.976 between the two).
- The median per-site difference is 0. Of the 4,996 sites (of 751,483) that differ by more than
  10 points at >= 10x, the median change in coverage between the BAMs is 26%, against 0% for the
  rest: they are sites where the two aligners placed different reads (repeats), not differences in
  how the same reads were scored. The 729 sites differing by more than 25 points are listed in
  `evidence/discordant_sites_gt25.tsv`.
- bonefish's pileup from the original BAM is byte-identical to that from the FASTQ: the BAM-input
  path restores read orientation and keeps the modification tags exactly.

## Time and memory

| | wall | CPU time | peak RSS |
|---|---:|---:|---:|
| minimap2 + samtools sort | 7 min 57 s | 5,149 s | 19.8 GB |
| bonefish `-y -o` (sorted BAM) | 2 min 50 s (1 min 51 s with the batch fix below) | 1,222 s | 23-26 GB |
| modkit pileup (each BAM) | 16-18 s | | |

Memory with methylation tags: before this test bonefish sized its read batches by bases only, and
ONT's tags take more memory than the bases, so `-y` held several large batches at once. Counting
the tags fixed it (`memory_*.time`): SAM output with `-y` peaked at 21.8 GB in v0.1.0 and 14.7 GB
after the fix, the same as without `-y` (15.0 GB). With sorted-BAM output the sorter's batches add
the rest (23-25 GB; 19.2 GB with `--sort-mem 1G`); per-second RSS of one run is in
`memory_fixed_sorted_y_rss_per_second.txt` (seconds, RSS, anonymous, file-backed, in kB).

## Files

- `evidence/comparison_cov{5,10,20}.tsv`: the comparison tables.
- `evidence/discordant_sites_gt25.tsv`: sites differing by > 25 points at >= 10x.
- `evidence/bam_input_check.txt`: FASTQ-input vs BAM-input pileups.
- `evidence/*.time`, `*.err`, `*.modkit.log`, `run.log`: timings, aligner and modkit logs.
- `evidence/memory_*`: the memory measurements above.
- `evidence/modkit_compare.sh`, `modkit_compare.py`: the scripts.
- `evidence/OUTPUTS_NOT_IN_GIT.sha256`: SHA-256 of the reads, BAMs and bedMethyl files kept outside Git.
- `SHA256SUMS` covers `evidence/`. Local paths are replaced by `${REPO}`, `${DATA}`, `${HOME}`, `${SCRATCH}`.
