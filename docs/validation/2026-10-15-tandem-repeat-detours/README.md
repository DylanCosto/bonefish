# Tandem-repeat false positives: cause and fix (15 Oct 2026)

The [tandem-repeat benchmark](../2026-10-14-cmrg-tandem-repeats/README.md) showed slightly more false
SV calls inside variable-number tandem repeats (VNTRs) from bonefish v0.1.3's BAMs than from minimap2's
(refined F1 0.9000 vs 0.9026 genome-wide). This traces the cause, and the fix in v0.1.4.

## Cause

At the false-positive loci, bonefish aligned reads carrying a small expansion of the repeat as a large
deletion plus a large insertion a few dozen bases apart, where minimap2 aligned a single small
insertion (`evidence/locus_examples.txt`):

| locus (chr20) | true allele | minimap2 | bonefish v0.1.3 |
|---|---|---|---|
| 53,500,100-53,500,450 | ~10 bp insertion | I10 in 19 of 40 reads | **D60 + I74** in 16 reads, I10 in 5 |
| 61,984,400-61,984,750 | ~13 bp insertion | I13 in 30 of 54 reads | **I76 + D63** in 23 reads, I13 in 5 |

Sniffles2 then reported the 60-76 bp gaps as separate SVs, both false. Switching off the parts of the
`dna` preset one at a time showed the source: **sub-anchored gap filling**. It splits long gaps at exact
15-mer matches that are unique in both the read and the reference window, and aligns only the pieces
between them. In a tandem repeat a match can be unique and still pair the read with another copy of
the repeat unit, so the alignment must detour through it: delete one copy's worth, then insert it
back. No other part of the preset changed these alignments.

Counting reads with such a detour (an insertion and a deletion of >= 30 bp within 300 bp, inside GIAB
Tier 1 tandem repeats; `evidence/detours.py`):

| | chr20 (55x) reads / loci | chr22 (27x) reads / loci |
|---|---:|---:|
| minimap2 2.31 | 74 / 21 | 10 / 10 |
| bonefish v0.1.3 | 405 / 79 | 138 / 49 |
| bonefish, sub-anchoring off | 76 / 23 | 10 / 10 |

## Options tested

HG002 chr20 (PAW70337, 55x) and chr22 (PAW71238, 27x), R10.4.1 sup; Clair3 + rtg vcfeval (GIAB
v4.2.1), Sniffles2 + Truvari (GIAB v5.0q), and the GIAB tandem-repeat benchmark v1.0.1 on each
chromosome's Tier 1 regions (Clair3 + Sniffles2 calls merged, `truvari bench` then `refine`):

| | TR F1, refined (chr20 / chr22) | TR false positives (sum) | SNP F1 (chr20 / chr22) | indel F1 | SVs vs v5.0q, TP / FP (sum) |
|---|---|---:|---|---|---|
| minimap2 2.31 | 0.8993 / 0.9030 | 512 | 0.9987 / 0.9963 | 0.9005 / 0.8854 | 960 / 100 |
| bonefish v0.1.3 | 0.8976 / 0.9008 | 531 | 0.9987 / 0.9963 | 0.9013 / 0.8853 | 961 / 96 |
| **sub-anchoring off (v0.1.4)** | 0.8983 / **0.9037** | 514 | 0.9987 / 0.9963 | 0.9010 / 0.8853 | 958 / 98 |
| sub-anchoring with detour repair | 0.8989 / 0.9014 | 516 | 0.9987 / 0.9964 | 0.9006 / 0.8854 | 957 / 101 |

The detour repair (experimental, `evidence/experiments.patch`) kept sub-anchoring but dropped the
segments around any detour and re-aligned that stretch. It removed the detours at no CPU cost, but
left a slightly worse SV benchmark on chr22 (61 false positives against 54).

**Decision: sub-anchoring is off in the `dna` preset from v0.1.4.** Alignments in tandem repeats are
then as clean as minimap2's and the tandem-repeat false positives fall to minimap2's level, with SNP
and indel accuracy unchanged and SV calls within 1-3 of v0.1.3 per chromosome. The option remains as
`--gap-subanchor`, which reproduces v0.1.3's output exactly.

## Costs and other checks

- **Speed:** on the genome-wide HG001 sample (16 threads, idle machine) 47.9 s and 698 CPU-s against
  44.7 s and 660 CPU-s with sub-anchoring: 5.8% more CPU, 7% more wall time. Measured again twice on
  a busy machine: 742-744 vs 707-708 CPU-s (+5%; `evidence/timing_busy_machine.txt`).
- **Placement:** on chr20 and chr22, no read becomes unmapped; 97.0% and 97.6% of primary alignments
  are identical and 0.14-0.15% move (`evidence/primary_changes.txt`). On the whole-genome simulated
  reads (`evidence/simulated/`), unmapped and wrongly placed counts are identical to v0.1.3 at every
  MAPQ threshold for all four read sets.
- **Output:** the genome-wide HG001 sample's MD5 changes from `deed7949d35f15bd8eed935be47a1d81`
  (v0.1.1-v0.1.3) to `8fbcd8ce9df67b78a978ef6497873caa`. The regression test data has no gaps long
  enough to be sub-anchored, so its expected output is unchanged.
- **Not re-run:** the whole-genome GIAB comparisons (their reads are no longer on disk). On chr20 and
  chr22 the change removes 17 of the 19 extra tandem-repeat false positives, which matches the
  genome-wide gap in size (+579 false positives over 24 chromosomes).

## Files

`evidence/score_chr20.tsv`, `score_chr22.tsv` (tandem-repeat benchmark and detour counts per run),
`runs/<run>/` (small-variant and SV accuracy, Truvari summaries, alignment time), `locus_examples.txt`,
`primary_changes.txt`, `timing_busy_machine.txt`, `simulated/`, the scripts (`trfp_runs.sh`, `trfp_runs2.sh`, `trfp_score.sh`,
`detours.py`) and `experiments.patch` (the diagnostic switch, the diagonal filter and the detour repair).
`SHA256SUMS` covers `evidence/`. Local paths are replaced by `${REPO}`, `${DATA}` and `${HOME}`.
