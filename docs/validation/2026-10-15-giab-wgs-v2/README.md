# Whole-genome re-validation after the tandem-repeat fix (15 Oct 2026)

The [12 Oct whole-genome validation](../2026-10-12-giab-wgs/README.md) was made with v0.1.1. Since
then, sub-anchored gap filling was found to misalign reads in tandem repeats
([v0.1.4 fix](../2026-10-15-tandem-repeat-detours/README.md)) and was brought back outside them
([this build](../2026-10-15-subanchor-outside-repeats/README.md)). Those changes were checked on
chr20 and chr22. This repeats the whole-genome runs with the new build.

## Setup

- **Build:** branch `subanchor-tr` at commit `e70c24e` (`dna` preset with `--gap-subanchor-tr 150`;
  the binary still reports version 0.1.4).
- **Reads:** the same as on 12 Oct, downloaded again: HG002 is 60% of flow cell PAW71238 (27x; same
  4,751,476 reads), HG005 all of PAW87816 (28x). R10.4.1, dorado sup v5.0.0.
- **Pipeline:** unchanged (`evidence/sample_eval.sh`): `bonefish -o aln.bam` (sorted, indexed),
  Clair3 2.0.3, rtg vcfeval against GIAB v4.2.1, Sniffles2 and Truvari against GIAB v5.0q SVs; for
  HG002 also v5.0q small variants (`q100_eval.sh`), CMRG v1.00 and the GIAB tandem-repeat benchmark
  v1.0.1 (`cmrg_tr_eval.sh`).
- **minimap2:** not re-run. The 12 Oct minimap2 2.31 alignments' calls and per-read tables were kept
  and are compared against as before (`wgs_v2.sh` links them).

## Results

| | HG005 minimap2 | HG005 v0.1.1 | HG005 **new** | HG002 minimap2 | HG002 v0.1.1 | HG002 **new** |
|---|---:|---:|---:|---:|---:|---:|
| SNP F1 (FP) | 0.9975 (8,774) | 0.9974 (8,915) | **0.9975** (8,887) | 0.9971 (11,752) | 0.9971 (11,851) | 0.9971 (11,785) |
| indel F1 | 0.9395 | 0.9396 | 0.9395 | 0.8737 | 0.8737 | 0.8737 |
| segmental-duplication SNP FP | 4,594 | 4,750 | 4,725 | 7,067 | 7,227 | 7,171 |
| SV F1 vs v5.0q, strict / multi-match | – | – | – | 0.8055 / 0.8706 | 0.8094 / 0.8741 | 0.8053 / 0.8706 |
| tandem repeats, refined F1 (FP) | – | – | – | 0.9026 (12,686) | 0.9000 (13,265) | **0.9023** (12,714) |
| CMRG SNP / indel / SV F1 | – | – | – | 0.9650 / 0.7514 / 0.9082 | 0.9649 / 0.7518 / 0.9032 | 0.9650 / 0.7518 / 0.9070 |
| primary alignments identical to minimap2 | | 90.6% | 90.6% | | 94.9% | 95.0% |
| reads newly unmapped vs minimap2 | | 3 | 3 | | 2 | 2 |
| reads -> sorted BAM, wall (CPU-s) | | 34.8 min (26,033) | 56.4 min* (25,013) | | 37.5 min (25,795) | **38.4 min** (24,351) |

\* HG005 was aligned while the HG002 download wrote to the same NTFS (ntfs-3g) drive: the sort's
merge stalled on writes (CPU ~17% for 20 minutes). HG002 ran alone and is the comparable timing.

HG002 against GIAB v5.0q small variants: SNP F1 0.9863 (minimap2 0.9862, v0.1.1 0.9863), indel F1
0.7075 (0.7075, 0.7076); outside difficult regions SNP 0.9972 and indel 0.9586, as v0.1.1.

## Structural variants and tandem repeats

The tandem-repeat gap of v0.1.1 is closed genome-wide: refined F1 0.9023 against minimap2's 0.9026,
with 28 more false positives than minimap2 instead of 579. For variants of 50 bp or more in Tier 1
repeats, false positives fall from 2,546 to 2,080 (minimap2 2,083) and missed variants from 2,394 to
2,335 (minimap2 2,302).

The SV benchmark against GIAB v5.0q moves the other way: strict F1 0.8094 -> 0.8053, now equal to
minimap2 (0.8055); 293 fewer true positives and 174 fewer false positives. Of the 717 truth SVs
v0.1.1 found and this build misses, 597 (83%) lie in Tier 1 tandem repeats; of the 424 found only by
this build, 320 do (`evidence/HG002/sv_lost_gained_vs_v011.txt`). Outside tandem repeats the
balance is 120 lost against 104 gained. Re-scoring with sequence comparison on (`--pctseq 0.7`
instead of 0) leaves all three results unchanged within 4 calls (`evidence/HG002/sv_sequence_matching/`),
but inside a repeat any copy of the unit resembles any other, so sequence similarity does not tell a
correct call from a detour-made one there. The tandem-repeat benchmark, which realigns the
haplotypes of disputed regions, finds more true large variants in repeats with this build than with
v0.1.1 (2,335 missed against 2,394). The v5.0q difference is therefore read as matching artefacts of
detour calls in repeats. Either way, bonefish's SV accuracy is now at minimap2's level rather than
above it.

## Files

Per sample (`evidence/HG002`, `evidence/HG005`): `report.md`, `run.log`, `resources.csv`,
`compare_primary.tsv`, `sv_agreement.tsv`, and in `bonefish/`: `accuracy.tsv`, `strat.tsv`,
`bam_summary.tsv`, `clair3_counts.tsv`, `flagstat.txt`, GNU time of each stage, rtg and Truvari
summaries, and for HG002 `q100.tsv`. `evidence/HG002/cmrg_tr/` (CMRG and tandem-repeat summaries),
`sv_sequence_matching/` (Truvari with sequence comparison, all three call sets) and
`sv_lost_gained_vs_v011.txt`. Scripts: `wgs_fetch.sh` and `s3_parallel_cat.py` (stream ONT's BAM to
FASTQ, optionally with parallel ranged reads), `wgs_v2.sh`, `sample_eval.sh`, `q100_eval.sh`,
`cmrg_tr_eval.sh`, `sv_lost_gained.sh`. The minimap2 and v0.1.1 inputs are in
[2026-10-12-giab-wgs](../2026-10-12-giab-wgs/README.md) and
[2026-10-14-cmrg-tandem-repeats](../2026-10-14-cmrg-tandem-repeats/README.md). `SHA256SUMS` covers
`evidence/`. Local paths are replaced by `${REPO}`, `${DATA}` and `${HOME}`.
