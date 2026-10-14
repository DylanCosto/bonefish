# Benchmarks and validation

All figures compare bonefish with minimap2 2.31 (`-ax lr:hq`, minimap2's preset for accurate
Nanopore reads) on the same reads, reference (GRCh38 no-alt analysis set) and machine. Accuracy is
judged the way the reads are used: small variants called with Clair3 and structural variants
with Sniffles2, scored against the GIAB truth sets.

**Machine.** Intel Core i5-12600K (16 threads, AVX2), 31 GB RAM, NVMe SSD, Linux; 16 threads for
both aligners. Apple Silicon results are in [validation/](validation/).

**Version.** The x86 runs below were made with builds whose output is byte-identical to v0.1.1
(checked on the genome-wide HG001 sample, MD5 `deed7949d35f15bd8eed935be47a1d81`). Since v0.1.1
arm64 builds produce the same output as x86 (same MD5 on an Apple M4 Pro).

## Whole genomes, reads to sorted and indexed BAM

ONT's GIAB 2025.01 data, R10.4.1, dorado sup v5.0.0. HG005: all reads of flow cell PAW87816
(27.9x); HG002: 60% of PAW71238 (27.2x). minimap2 `-ax lr:hq | samtools sort -@16 --write-index`
against `bonefish -o sample.bam` (built-in sort and index). Clair3 2.0.3 (sup v5.0.0 model) and
rtg vcfeval against GIAB v4.2.1; Sniffles2 and Truvari against GIAB v5.0q (HG002). Per-sample
reports and raw results: [validation/2026-10-12-giab-wgs](validation/2026-10-12-giab-wgs/README.md).

| | HG005 minimap2 | HG005 bonefish | HG002 minimap2 | HG002 bonefish |
|---|---:|---:|---:|---:|
| reads to sorted, indexed BAM | 114.7 min | **34.8 min (3.3x)** | 114.8 min | **37.5 min (3.1x)** |
| CPU time | 25.3 h | 7.2 h | 25.3 h | 7.2 h |
| peak memory (process) | 16.7 GB + sort 7.3 GB | 19.6 GB | 16.9 GB + sort 7.2 GB | 19.7 GB |
| SNP F1 | 0.9975 | 0.9974 | 0.9971 | 0.9971 |
| indel F1 | 0.9395 | 0.9396 | 0.8737 | 0.8737 |
| SV F1 vs v5.0q, strict / multi-match | – | – | 0.8055 / 0.8706 | **0.8094 / 0.8741** |
| primary alignments identical / moved | | 90.6% / 1.1% | | 94.9% / 1.4% |
| reads newly unmapped | | 3 | | 2 |

HG002 small variants re-scored against GIAB v5.0q (the T2T-Q100-based truth, 2.74 Gbp including
chrX/Y and difficult regions):

| region | SNP F1 minimap2 (FP) | SNP F1 bonefish (FP) | indel F1 minimap2 | indel F1 bonefish |
|---|---:|---:|---:|---:|
| all | 0.9862 (46,658) | 0.9863 (45,136) | 0.7075 | 0.7076 |
| new vs v4.2.1 | 0.8774 (42,180) | 0.8787 (40,782) | 0.4475 | 0.4475 |
| segmental duplications | 0.9333 (8,857) | 0.9337 (8,751) | 0.6506 | 0.6512 |
| not in difficult regions | 0.9971 | 0.9972 | 0.9584 | 0.9586 |

## Methylation calls

HG002 chr20 (55x, dorado sup v5.0.0 with 5mC/5hmC models), aligned by both with `-y`, then
`modkit pileup --cpg --combine-strands`: per-CpG 5mC agrees with minimap2's at r = 0.9988 (sites
>= 10x in both), mean absolute difference 0.21 percentage points, chr20 5mC level 64.42% vs 64.38%;
the few discordant sites are where the aligners placed different reads. Details:
[validation/2026-10-14-modkit-methylation](validation/2026-10-14-modkit-methylation/README.md).

## Aligner speed alone (SAM output)

| dataset | minimap2 | bonefish | speedup |
|---|---:|---:|---:|
| HG001 R10.4.1 sup, random genome-wide sample (2.7 Gbp, 161,762 reads) | 202 s | 44 s | 4.6x |
| HG002 R10.4.1 HAC, chr20 reads (3.5 Gbp) | 417 s | 72 s | 5.8x |
| HG002 R9.4.1, chr20 reads (4.6 Gbp) | 412 s | 103 s | 4.0x |
| HG002 R10.4.1 ultra-long (3.4 Gbp) | 542 s | 65 s | 8.3x |
| *Salmonella enterica* R10.4.1 (1.2 Gbp) | 25 s | 10 s | 2.6x |

## Known differences and costs

- Output is not byte-identical to minimap2: 1–2% of primary alignments are placed differently.
- Segmental duplications against GIAB v4.2.1: about 160 more SNP false positives per genome
  (F1 −0.0007); against v5.0q bonefish has fewer false positives there.
- About 0.4% fewer read bases are aligned: low-identity read ends are left soft-clipped (primary
  alignment identity is higher, 0.979 vs 0.977).
- Clair3 takes 9–12% longer on bonefish BAMs.
- Peak memory for a human genome is about 20 GB with 16 threads (the GRCh38 index is 7.3 GB, or
  8.5 GB in the raw format that loads in under a second).

## Other ONT aligners (14 Oct 2026)

Raw timings, logs, statistics and the script: [validation/2026-10-14-aligner-comparison](validation/2026-10-14-aligner-comparison/README.md).

Same machine, reference, reads and 16 threads; one run each; SAM output unless noted. Each tool
uses its own index and its recommended ONT preset: `lr:hq` for minimap2, dorado (default) and
mm2plus; `map-ont` for mm2-fast (based on minimap2 2.24, which predates `lr:hq`) and Winnowmap2;
`-x ont` for NGMLR, `-ONT` for LRA. mm2-fast ran its AVX2 code (this CPU has no AVX-512).

**Genome-wide HG001 sample** (161,762 R10.4.1 sup reads, 2.7 Gbp):

| aligner | wall | vs bonefish | CPU time | peak memory | reads mapped | bases aligned | identity |
|---|---:|---:|---:|---:|---:|---:|---:|
| **bonefish 0.1.1** | **46 s** | 1.0x | 658 s | 15.0 GB | 96.30% | 96.44% | 0.979 |
| minimap2 2.31 | 192 s | 4.2x | 2,563 s | 18.3 GB | 96.30% | 96.78% | 0.977 |
| dorado 2.1.2 aligner | 198 s | 4.3x | 2,167 s | 13.3 GB | 96.30% | 95.85% | 0.979 |
| mm2plus 1.3 | 236 s | 5.1x | 3,059 s | 18.5 GB | 96.30% | 96.78% | 0.977 |
| mm2-fast (minimap2 2.24) | 641 s | 13.9x | 9,488 s | 19.6 GB | 96.99% | 97.14% | 0.974 |
| Winnowmap2 2.03 | 1,764 s | 38x | 25,356 s | 19.4 GB | 96.19% | 96.83% | 0.980 |

Winnowmap 2.03 maps nothing from a prebuilt index, so its time includes building one (63 s).

**10% of the same reads** (16,177 reads, 270 Mbp), to include the slower aligners (times include
loading the index, which favours tools with small or fast-loading indexes on short runs):

| aligner | wall | CPU time | reads mapped |
|---|---:|---:|---:|
| **bonefish** | **11 s** | 68 s | 96.30% |
| dorado aligner | 27 s | 210 s | 96.30% |
| mm2plus | 34 s | 299 s | 96.30% |
| minimap2 | 34 s | 268 s | 96.30% |
| mm2-fast | 68 s | 932 s | 96.92% |
| LRA 1.3.7 | 75 s | 702 s | 89.56% |
| Winnowmap2 | 203 s | 2,307 s | 96.08% |
| NGMLR 0.2.7 | 2,909 s | 24,561 s | 93.00% |

**Reads to sorted, indexed BAM** (genome-wide sample; output on the same disk):

| pipeline | wall | CPU time | peak memory |
|---|---:|---:|---:|
| **bonefish `-o out.bam`** | **55 s** | 758 s | 16.4 GB |
| minimap2 `-ax lr:hq` \| samtools sort `--write-index` | 219 s (4.0x) | 2,699 s | 18.6 GB |
| dorado aligner `-o DIR` (sorted, indexed) | 574 s (10.5x) | 2,501 s | 13.3 GB |

**Index build:** bonefish 62 s; minimap2, mm2plus and mm2-fast 61–65 s; Winnowmap 63 s plus
130 s of k-mer counting (meryl); LRA 196 s; NGMLR 416 s.

Notes: the mapping columns are not accuracy (see the variant-calling results above for
bonefish against minimap2). dorado aligner reports far more supplementary alignments than minimap2
on the same index (22,430 vs 13,543), so its defaults differ from minimap2's command line.
Identity is over aligned bases, so tools that clip more score higher; LRA and NGMLR leave 7-10% of
reads unmapped. Swap use rose briefly during the minimap2 run (to 7.9 GB); its time matches earlier
measurements on this machine.
