# Sub-anchoring outside tandem repeats (15 Oct 2026)

v0.1.4 turned off sub-anchored gap filling, because inside tandem repeats it aligned reads through
the wrong copy of the repeat unit, which callers reported as extra SVs
([details](../2026-10-15-tandem-repeat-detours/README.md)). That cost about 6% CPU. This change
turns sub-anchoring back on everywhere except in tandem repeats.

## Change

Sub-anchoring splits a long gap at exact 15-mer matches that are unique in both the read and the
reference window, and aligns only the pieces between them. It now drops every such segment that
overlaps a tandem repeat in the read or the reference. A tandem repeat is marked wherever a 15-mer
recurs within 150 bp; the mark spans from its first copy to the end of its last. The stretch
around the repeat is then aligned in one piece, as v0.1.4 aligns the whole gap. The `dna` preset
sets `--gap-subanchor-tr 150`. `--gap-subanchor-tr 0` gives v0.1.3's output and
`--no-gap-subanchor` gives v0.1.4's, both byte for byte.

## Choosing the repeat period

Reads with a deletion-insertion detour (an insertion and a deletion of >= 30 bp within 300 bp)
inside GIAB Tier 1 tandem repeats, HG002 (`evidence/detour_screen.tsv`; the screen used an
environment variable that became `--gap-subanchor-tr`, with identical output):

| largest period masked | chr22 detour reads | chr22 CPU-s | chr20 detour reads | chr20 CPU-s |
|---|---:|---:|---:|---:|
| v0.1.3 (no mask) | 138 | 608 | 405 | 1,081 |
| 40 bp | 80 | 608 | | |
| 80 bp | 14 | 616 | | |
| **150 bp** | **10** | 616 | **76** | 1,096 |
| 200 bp | | | 75 | 1,126 |
| 300 bp | 10 | 623 | | |
| v0.1.4 (sub-anchoring off) | 10 | 648 | 76 | 1,175 |
| minimap2 2.31 | 10 | | 74 | |

150 bp is the shortest period that removes every extra detour on both chromosomes.

## Speed

Genome-wide HG001 sample, 16 threads, two interleaved runs each (`evidence/timing_genome_wide_sample.txt`;
a download used about half a core throughout):

| | wall | CPU-s |
|---|---:|---:|
| v0.1.3 | 46.6 s | 664 |
| **this change** | **46.9 s** | **671** |
| v0.1.4 | 49.7 s | 700 |

80% of v0.1.4's extra CPU time is recovered; the remaining cost against v0.1.3 is 1%.

## Accuracy

HG002 chr20 (PAW70337, 55x) and chr22 (PAW71238, 27x), R10.4.1 sup, as for v0.1.4:

| | TR F1, refined (chr20 / chr22) | detour reads | SNP F1 | indel F1 | SVs vs v5.0q, TP / FP (chr20; chr22) |
|---|---|---|---|---|---|
| minimap2 2.31 | 0.8993 / 0.9030 | 74 / 10 | 0.9987 / 0.9963 | 0.9005 / 0.8854 | |
| v0.1.3 | 0.8976 / 0.9008 | 405 / 138 | 0.9987 / 0.9963 | 0.9013 / 0.8853 | 566 / 44; 395 / 52 |
| v0.1.4 | 0.8983 / 0.9037 | 76 / 10 | 0.9987 / 0.9963 | 0.9010 / 0.8853 | 564 / 44; 394 / 54 |
| **this change** | **0.8984 / 0.9046** | 76 / 10 | 0.9987 / **0.9965** | 0.8992 / 0.8845 | 564 / 45; 394 / **44** |

- **Tandem repeats and SVs:** as v0.1.4 on chr20. On chr22 the tandem-repeat F1 is the highest of
  the four and SV false positives the fewest (10 fewer than v0.1.4).
- **Indels:** F1 is 0.0018 / 0.0008 lower than v0.1.4. The changed calls are Clair3 churn, not an
  alignment effect: of the 72 new false positives and the 154 indels v0.1.3 called that this run
  misses, 15% and 14% have any primary-alignment change within 300 bp, against 16% of randomly
  chosen truth indels (`evidence/indel_churn_chr20.tsv`). Clair3 phases the whole chromosome, so
  alignment changes anywhere shift calls elsewhere. Across about 60 earlier chr20 runs of alignment
  variants, indel F1 ranged from 0.8996 to 0.9016.
- **Placement:** against v0.1.4, 99.2% (chr20) and 99.6% (chr22) of primary alignments are
  identical and 0.04-0.08% move; no read becomes unmapped (`evidence/primary_changes.txt`).
- **Simulated reads:** whole-genome simulations (`evidence/simulated/`): unmapped and wrongly placed
  counts identical to v0.1.4 at every MAPQ threshold for all four read sets.

## Files

`evidence/runs/<run>/` (small-variant and SV accuracy, Truvari tandem-repeat summaries, alignment
time), `score_chr20.tsv` and `score_chr22.tsv` (tandem-repeat benchmark and detours, all runs),
`detour_screen.tsv`, `timing_genome_wide_sample.txt`, `indel_churn_chr20.tsv`, `primary_changes.txt`,
`simulated/`, and the scripts (`subtr_quick.sh`, `subtr_runs.sh`, `indel_churn.sh`, `detours.py`,
`trfp_score.sh`). `SHA256SUMS` covers `evidence/`. Local paths are replaced by `${REPO}`, `${DATA}`
and `${HOME}`.
