# Spliced alignment: direct RNA and cDNA (21 Oct 2026)

bonefish was DNA-only from 11 Oct. This adds the spliced presets back and speeds them up:

- `-x rna`: direct RNA, minimap2 2.31's `-ax splice -uf -k14`.
- `-x cdna`: cDNA, minimap2's `-ax splice`. `splice` is an alias, and `splice:hq` is minimap2's
  accurate-cDNA preset.

The restored code reproduces minimap2 byte for byte. The speedups change about 0.1-0.5% of
primary alignments. The `dna` preset's output is unchanged.

## Data

- **HG002 direct RNA:** RNA004 chemistry, dorado sup, 902,275 reads (1.14 Gbp, mean 1.3 kb),
  with m6A and other modification tags. **HG002 cDNA:** PCS114 kit, hac, 808,749 reads (2.04 Gbp,
  mean 2.5 kb). Both are the first 2 GB of the unaligned BAMs of ONT's `UHRR_HG002_2026.06`
  release (`evidence/scripts/rna_fetch.sh`). Iterations used the first 200,000 reads of each.
- **Reference:** GRCh38 no-alt analysis set. **Annotation:** GENCODE v49 (as in ONT's own analysis
  of these runs), 507,365 transcripts as BED12.
- **Simulated reads with known intron chains:** 30,000 random multi-exon GENCODE v49 transcripts
  (`evidence/scripts/rna_sim_reads.sh`). Badread reads like cDNA (98% mean identity, both
  strands, 44,667 reads) and like direct RNA (95%, sense strand, 44,787 reads). A read counts as
  correct ("exact intron chain") when its primary alignment has exactly the true introns
  (`rna_sim_truth.py`).
- **Comparison:** minimap2 2.31 with the same presets. Machine: Core i5-12600K, 16 threads, 31 GB.
  Both aligners used prebuilt indexes, paged in before each run.

## Speed

| | minimap2 | bonefish (restored) | **bonefish** |
|---|---:|---:|---:|
| direct RNA, 200k reads: CPU-s | 1,215 | 1,164 | **564** |
| cDNA, 200k reads: CPU-s | 4,298 | 4,003 | **795** |
| direct RNA, all 902k reads: CPU-s | 5,606 | | **2,618** |
| cDNA, all 809k reads: CPU-s | 17,416 | | **3,275** |
| direct RNA, all reads: wall, idle machine | 455.6 s | | **177.3 s** |
| cDNA, all reads: wall, idle machine | 1,222.8 s | | **227.0 s** |

- **Speed-up over minimap2:** 2.1x less CPU time on direct RNA and 5.3x on cDNA. bonefish
  keeps more cores busy, so wall time falls more: 2.6x and 5.4x on the full sets.
- **Wall-time caveat:** the full-set wall times come from the build before the last commit's
  rescue fix, which costs +1% / +3% CPU. Later runs shared the machine with interactive use,
  so their wall times are not comparable.

## Accuracy

| | minimap2 | bonefish |
|---|---:|---:|
| direct RNA (all reads): primaries identical to minimap2 | | 99.87% |
| direct RNA: same intron chain as minimap2 (reads spliced by either) | | 99.78% |
| direct RNA: mapped | 80.70% | 80.70% |
| direct RNA: introns at GENCODE-annotated positions | 97.94% | 97.95% |
| cDNA (all reads): primaries identical to minimap2 | | 99.49% |
| cDNA: same intron chain as minimap2 | | 99.44% |
| cDNA: mapped | 99.37% | 99.37% |
| cDNA: introns at GENCODE-annotated positions | 97.49% | 97.46% |
| simulated direct-RNA-like: exact intron chain | 89.10% | 88.79% |
| simulated direct-RNA-like: intron precision / recall | 98.48% / 97.34% | 98.37% / 97.09% |
| simulated cDNA-like: exact intron chain | 93.13% | 93.09% |
| simulated cDNA-like: intron precision / recall | 97.57% / 97.81% | 97.56% / 97.80% |

On the cDNA set, the transcript-strand calls (`ts:A`) match minimap2's to within 0.05%:
`+` 90,650 vs 90,635, `-` 91,905 vs 91,871.

## What changed

1. **Intron gaps.** Profiling the restored code showed where the time went. Gap fills between
   chain anchors that span an intron took 85-86% of spliced alignment time: about 2.3 million
   DP cells per intron, because the splice DP fills the whole query x reference matrix.
   - **Why it can be cut:** the intron's second gap cost has no per-base extension, so its cost
     is the same at any length.
   - **Compression:** reference stretches far from every chain anchor are replaced by a 256-base
     N seam. The query is aligned to the compressed reference, and each intron that crosses a
     seam is widened back to its real length.
   - **Safety check:** a path that crosses a seam any other way than inside one intron falls
     back to the full DP.
   - **Missing exons:** an exon with no anchor inside the removed reference can't be aligned this
     way. It shows up as poor alignment next to the widened intron: more than 6 errors within 32
     query bases. In that case, exact 10-mer matches of that query stretch are searched for in
     the removed reference, and windows around the best 6 candidates are kept for a second
     alignment (two rounds). Without this search, the simulated direct-RNA exact-intron-chain
     rate was 86.3% instead of 88.8%.
2. **Strands.** For cDNA (`-ub`), minimap2 aligns each chain twice, once per transcript strand.
   bonefish aligns it once, scoring every splice site by whichever strand's motif fits better.
   That score is an upper bound on both strands.
   - **One strand wins:** if all the path's introns score as well under one strand as under the
     bound, and worse under the other, the path is that strand's best alignment.
   - **Tie:** with no introns, or introns that score equally under both strands, the strands tie,
     as minimap2 reports.
   - **Fallback:** mixed cases, and every run with `--junc-bed`, align both strands as before.
3. **Chaining.** The spliced chainer now uses the SIMD chainers (AVX-512, AVX2, SSE and NEON).
   It has minimap2's intron gap cost and gives output identical to the scalar chainer.

## Tried and not kept

- **Full-DP fallback on poor flanks** instead of the exon search: accuracy matched minimap2, but
  it cost +25-45% CPU.
- **Lower seed occurrence cap** (`-M 1000` / `-M 500` for direct RNA): 16-27% faster, but 337 /
  477 of 200k reads became unmapped (repeats and paralogs).
- **Stable LSD radix sort of anchors:** 40% faster seeding on one thread, but +6% CPU at 16
  threads, where sorting about 21,000 anchors per read is limited by memory bandwidth.
- **Dropping isolated anchors before chaining:** output identical, no measurable change.

Details and intermediate measurements: `evidence/optimisation_log.md`, `evidence/ledger.tsv`.

## Known limits

Both are documented in the tests.

- **Exon with no exact 10-mer:** when an exon has no anchor and none of its 10-mers matches the
  reference exactly, the compressed alignment misses it and keeps a worse alignment. This hit
  1 of 225 cases in a random-gene stress test at about 8% read error.
- **Hidden ties:** the single cDNA pass cannot see a tie between strands reached through a
  different path (1 of 400 random reads), so it reports a strand where two passes would report
  none.

## DNA

The `dna` preset is unchanged:

- **Golden-output regression tests** pass without re-blessing.
- **Genome-wide HG001 sample:** identical SAM to v0.1.7 (MD5 `9bc7ff3e9f498cc1530b2b1d78c47e10`).
- **Speed:** 650-651 vs 648 CPU-s over three alternating rounds (+0.4%, `evidence/dna_guard.txt`).

## Files

- `evidence/runs/<run>/`: GNU time, stage timers, comparisons with minimap2
  (`compare_primary.py`, `intron_chains.py`) and simulated-read scores, for each run:
  - `mm2`, `mm2full`: minimap2;
  - `base`: the restored code;
  - `c8_s`, `optfull`: the build before the rescue fix;
  - `final`, `final_full`: this branch.
- `evidence/ledger.tsv`: one line per run.
- `evidence/sim_transcripts_truth.tsv.gz`: exon blocks of the simulated transcripts.
- `evidence/scripts/`: `rna_fetch.sh`, `rna_sim_reads.sh`, `rna_sim_transcripts.py`,
  `rna_sim_truth.py`, `rna_bench.sh`, `intron_chains.py`, `compare_primary.py`, `dna_guard.sh`.

`SHA256SUMS` covers `evidence/`. Local paths are replaced by `${BENCH}`, `${DATA}`, `${REPO}` and
`${HOME}`.
