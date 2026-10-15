# Medically relevant genes and tandem repeats, HG002 whole genome (14 Oct 2026)

Two GIAB benchmarks that target the regions short-read and older long-read pipelines get wrong:

- **CMRG v1.00** (Challenging Medically Relevant Genes): 273 genes of clinical interest that are hard
  to map (segmental duplications, gene families, repeats), with small-variant and SV truth from a
  curated diploid assembly.
- **GIAB tandem-repeat benchmark v1.0.1** (English et al.): variants of 5 bp or more in 1.7 million
  tandem-repeat regions; scored on the 1.64 million Tier 1 regions, as its README recommends.

## Data and setup

- **Calls:** the whole-genome HG002 calls of the [GIAB whole-genome comparison](../2026-10-12-giab-wgs/README.md)
  (60% of flow cell PAW71238, 27x, R10.4.1 sup; Clair3 2.0.3 and Sniffles2 on the minimap2 2.31 and
  bonefish BAMs). Nothing was realigned or re-called.
- **CMRG small variants:** rtg vcfeval of the PASS SNPs and indels against `HG002_GRCh38_CMRG_smallvar_v1.00`
  within its BED. **CMRG SVs:** Truvari 5.5 against `HG002_GRCh38_CMRG_SV_v1.00` within its BED, with
  `--pctseq 0 --typeignore --dup-to-ins` (the CMRG README's advice: its duplications are encoded as
  insertions).
- **Tandem repeats:** the Clair3 PASS calls and Sniffles2 PASS SVs merged into one call set
  (multi-allelic sites split), then `truvari bench --sizemin 5 --pick ac` on the Tier 1 regions and
  `truvari refine --use-original-vcfs`, which realigns the haplotypes of disputed regions with MAFFT so
  differently represented but equivalent calls count as matches (the benchmark's recommended method).
- Script: `evidence/cmrg_tr_eval.sh`.

## Results

| benchmark | minimap2 TP / FP / FN | minimap2 F1 | bonefish TP / FP / FN | bonefish F1 |
|---|---|---:|---|---:|
| CMRG SNPs | 17,334 / 964 / 295 | 0.9650 | 17,304 / 933 / 325 | 0.9649 |
| CMRG indels | 2,721 / 895 / 907 | 0.7514 | 2,722 / 893 / 906 | **0.7518** |
| CMRG SVs | 193 / 16 / 23 | **0.9082** | 196 / 22 / 20 | 0.9032 |
| tandem repeats, Tier 1 (bench) | 114,134 / 21,639 / 19,303 | **0.8469** | 113,883 / 22,473 / 19,554 | 0.8432 |
| tandem repeats, Tier 1 (refined) | 120,081 / 12,686 / 13,356 | **0.9026** | 119,895 / 13,265 / 13,542 | 0.9000 |

- **Medically relevant genes:** the same accuracy. bonefish has fewer false-positive SNPs (933 vs 964)
  and 30 more missed ones; it finds 3 more true SVs and makes 6 more false SV calls.
- **Tandem repeats:** bonefish is slightly behind (refined F1 0.9000 vs 0.9026). By size, in the bench
  output:

  | variants in Tier 1 repeats | minimap2 FP / FN | bonefish FP / FN |
  |---|---|---|
  | 5-49 bp (mostly Clair3 calls) | 19,536 / 17,001 | 19,904 / 17,160 |
  | 50 bp and over (mostly Sniffles2 calls) | 2,083 / 2,302 | 2,546 / 2,394 |

  Most of the difference is larger variants in VNTRs (+463 false positives; 16,149 truth
  variants of that size), where read placement and alignment through long repeats decide how the
  caller sees an allele. Genome-wide, SV accuracy against GIAB v5.0q is slightly higher from
  bonefish's BAM (F1 0.8094 vs 0.8055), so this is specific to tandem repeats. It is a known
  difference, listed in the [benchmarks](../../benchmarks.md#known-differences-and-costs).

  **Follow-up (15 Oct):** traced to sub-anchored gap filling aligning reads through the wrong copy of
  a repeat, and fixed in v0.1.4 ([details](../2026-10-15-tandem-repeat-detours/README.md)).

## Files

- `evidence/summary.tsv`: the table above.
- `evidence/<aligner>/`: rtg vcfeval and Truvari summaries, Truvari parameters, and
  `tr_refine_regions.txt.gz` (per-region results of the refine step).
- `evidence/cmrg_tr_eval.sh`; `TRUTH_FILES.sha256` (the GIAB files); `CALLS_NOT_IN_GIT.sha256` (the
  call sets scored).
- `SHA256SUMS` covers `evidence/`. Local paths are replaced by `${REPO}`, `${DATA}` and `${HOME}`.
