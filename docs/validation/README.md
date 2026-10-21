# Validation reports

Each folder holds the README for one validation round, its evidence (logs, tables, scripts) and a
`SHA256SUMS` file. [benchmarks.md](../benchmarks.md) summarises them.

## Accuracy against truth sets (DNA)

| report | what |
|---|---|
| [2026-10-12-giab-wgs](2026-10-12-giab-wgs/README.md) | Whole genomes HG002 and HG005: Clair3, Sniffles2, GIAB v4.2.1 and v5.0q (v0.1.1) |
| [2026-10-15-giab-wgs-v2](2026-10-15-giab-wgs-v2/README.md) | The same, re-run with v0.1.5 |
| [2026-10-14-hac-r9-phasing](2026-10-14-hac-r9-phasing/README.md) | HAC and R9.4.1 reads, phasing |
| [2026-10-14-cmrg-tandem-repeats](2026-10-14-cmrg-tandem-repeats/README.md) | 273 challenging medically relevant genes, GIAB tandem-repeat benchmark |
| [2026-10-15-tandem-repeat-detours](2026-10-15-tandem-repeat-detours/README.md) | Tandem-repeat SV false positives traced and fixed (v0.1.4) |
| [2026-10-15-subanchor-outside-repeats](2026-10-15-subanchor-outside-repeats/README.md) | Sub-anchoring restored outside tandem repeats (v0.1.5) |
| [2026-10-14-simulated-reads](2026-10-14-simulated-reads/README.md) | Simulated reads with known origin (R10-like, ultra-long, 0.5 kb, 93% identity) |
| [2026-10-14-modkit-methylation](2026-10-14-modkit-methylation/README.md) | CpG methylation calls with modkit |
| [2026-10-14-wf-human-variation](2026-10-14-wf-human-variation/README.md) | ONT's wf-human-variation pipeline |
| [2026-10-14-other-references](2026-10-14-other-references/README.md) | T2T-CHM13 and GRCh38 with ALT contigs |

## Other read types

| report | what |
|---|---|
| [2026-10-21-rna-cdna](2026-10-21-rna-cdna/README.md) | Direct RNA and cDNA (spliced presets) against minimap2, GENCODE and simulated reads |
| [2026-10-21-transcriptome](2026-10-21-transcriptome/README.md) | Transcriptome mapping for quantification (`-x tx`, oarfish) |
| [2026-10-20-metagenomics-duplex-rdna](2026-10-20-metagenomics-duplex-rdna/README.md) | Metagenome (Zymo fecal), HG002 duplex reads, rDNA |

## Speed, platforms and output

| report | what |
|---|---|
| [2026-10-14-aligner-comparison](2026-10-14-aligner-comparison/README.md) | Against dorado, mm2plus, mm2-fast, Winnowmap2, NGMLR and LRA |
| [2026-10-14-scaling-memory](2026-10-14-scaling-memory/README.md) | Threads and memory |
| [2026-10-17-avx512](2026-10-17-avx512/README.md) | AVX-512 kernels |
| [2026-10-13-macos-arm64](2026-10-13-macos-arm64/README.md), [2026-10-13-neon-clean-dp](2026-10-13-neon-clean-dp/README.md), [2026-10-13-hg001-macos-parity](2026-10-13-hg001-macos-parity/README.md) | Apple Silicon: NEON kernels, same output as x86 |
| [2026-10-14-output-checks](2026-10-14-output-checks/README.md) | SAM/BAM conformance (Picard), determinism across threads/inputs/indexes, large references |
| [2026-10-16-bam-index-chunks](2026-10-16-bam-index-chunks/README.md) | BAM index identical to samtools' (Clair3 speed fix, v0.1.6) |

Scripts to reproduce the platform checks are in [scripts/validation](../../scripts/validation/README.md).
