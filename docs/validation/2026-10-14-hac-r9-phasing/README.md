# HAC and R9.4.1 reads, and phasing (14 Oct 2026)

The earlier variant-calling checks used R10.4.1 sup reads. Many users basecall with the faster HAC
model, and much existing data is R9.4.1. This compares bonefish with minimap2 2.31 (`-ax lr:hq`) on
all three for small variants, structural variants and phasing.

## Data and setup

All reads are HG002 reads that map to chr20 (3.5-4.6 Gbp per set, 55-70x of chr20), aligned
against the whole GRCh38 no-alt genome:

| set | reads | basecaller / Clair3 model |
|---|---|---|
| R10.4.1 sup | flow cell PAW70337 (GIAB 2025.01) | dorado sup v5.0.0 / `r1041_e82_400bps_sup_v500` |
| R10.4.1 HAC | the same flow cell (3.47 Gbp) | dorado hac v5.0.0 / `r1041_e82_400bps_hac_v500` |
| R9.4.1 sup | ONT open data, HG002 R9.4.1 (4.56 Gbp) | Guppy 5.0.6 sup / `r941_prom_sup_g5014` |

- **Pipeline** (`evidence/eval_accuracy.sh`, unchanged from the earlier rounds): align and sort, Clair3
  2.0.3 on chr20, rtg vcfeval against GIAB v4.2.1 (PASS SNPs and indels, chr20 benchmark regions);
  Sniffles2 and Truvari against GIAB v5.0q (`evidence/eval_sv.sh`).
- **Phasing** (`evidence/phasing.sh`), on the sup BAMs: longphase (the phaser in ONT's
  wf-human-variation) phases the Clair3 PASS SNPs and indels; whatshap compare scores the phased SNVs
  against GIAB v5.0q, which comes from the fully phased Q100 assembly, within its benchmark regions;
  whatshap stats gives block lengths; longphase haplotag assigns reads to haplotypes.
- bonefish is this release's code (alignment output identical to v0.1.2); 16 threads.

## Results

| set | SNP F1 minimap2 | SNP F1 bonefish | indel F1 minimap2 | indel F1 bonefish | SV F1 minimap2 | SV F1 bonefish |
|---|---:|---:|---:|---:|---:|---:|
| R10.4.1 sup | 0.9987 | 0.9987 | 0.9005 | **0.9013** | 0.8011 | **0.8017** |
| R10.4.1 HAC | 0.9985 | 0.9984 | 0.8294 | **0.8310** | 0.7887 | **0.7917** |
| R9.4.1 sup | 0.9980 | 0.9980 | 0.8173 | 0.8172 | 0.7704 | **0.7734** |

True positives, false positives and false negatives are in `evidence/runs/<set>_<aligner>/`. SNP
true positives differ by at most 2 per set (of ~71,300); bonefish finds 2-8 more true SVs per set.
The minimap2 numbers reproduce the 11 Oct runs of the same pipeline exactly.

| phasing (sup, chr20) | minimap2 | bonefish |
|---|---:|---:|
| phased SNVs compared with the truth | 47,035 | 47,021 |
| switch errors (rate) | 7 (0.015%) | 7 (0.015%) |
| blockwise Hamming errors | 2,911 | 2,911 |
| phase blocks / block N50 | 64 / 2.17 Mb | 60 / 2.20 Mb |
| reads assigned to a haplotype | 76.0% | 76.2% |

Phasing is the same from either aligner's reads. The Hamming errors come from the few switches
inside long blocks: a switch in the middle of a block flips every variant after it, so both aligners
have the same count from the same switches.

**Speed of the alignment step** (reads to sorted BAM, `align.time`; Clair3 ran afterwards):

| set | minimap2 + samtools sort | bonefish + samtools sort |
|---|---:|---:|
| R10.4.1 sup | 369 s | 65 s (5.7x) |
| R10.4.1 HAC | 371 s | 65 s (5.7x) |
| R9.4.1 sup | 451 s | 109 s (4.1x) |

## Notes

- The simulated low-accuracy reads (~93% identity) in the
  [simulation check](../2026-10-14-simulated-reads/README.md) map better with minimap2's `map-ont`
  preset. Real R9.4.1 sup reads give the same accuracy from bonefish as from minimap2
  `lr:hq`, as the table shows.
- Clair3 takes 5-8% longer on bonefish's BAMs, as on whole genomes.
- One sample, one chromosome per set.

## Files

- `evidence/runs/<set>_<aligner>/`: `accuracy.tsv` (SNPs, indels), `sv_accuracy.tsv`, the rtg and
  Truvari summaries, `sv_calls.vcf.gz` (Sniffles2), `align.time`, `clair3.time`, `eval.log`.
- `evidence/phasing/`: `summary.tsv`, and per aligner the whatshap compare and stats tables and the
  longphase logs.
- `evidence/eval_accuracy.sh`, `eval_sv.sh`, `phasing.sh`, `validation_queue.sh` (the order the runs
  were made in), `queue.log`.
- `evidence/OUTPUTS_NOT_IN_GIT.sha256`: SHA-256 of the Clair3 VCFs.
- `SHA256SUMS` covers `evidence/`. Local paths are replaced by `${REPO}`, `${DATA}` and `${HOME}`.
