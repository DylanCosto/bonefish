# Benchmarks and validation

All figures compare bonefish with minimap2 2.31 (`-ax lr:hq`, minimap2's preset for accurate
Nanopore reads) on the same reads, reference (GRCh38 no-alt analysis set) and machine. Accuracy is
judged the way the reads are used: small variants called with Clair3 and structural variants
with Sniffles2, scored against the GIAB truth sets.

**Machine.** Intel Core i5-12600K (16 threads, AVX2), 31 GB RAM, NVMe SSD, Linux; 16 threads for
both aligners. Apple Silicon results are in [validation/](validation/).

**Version.** The x86 runs below were made with builds whose output is byte-identical to v0.1.1
(checked on the genome-wide HG001 sample, MD5 `deed7949d35f15bd8eed935be47a1d81`). Since v0.1.1
arm64 builds produce the same output as x86 (same MD5 on an Apple M4 Pro). v0.1.4 turned off
sub-anchored gap filling (about 6% more CPU) to fix the tandem-repeat false positives below
([details](validation/2026-10-15-tandem-repeat-detours/README.md)); v0.1.5 sub-anchors again outside
tandem repeats (1% more CPU than v0.1.3, MD5 `9bc7ff3e9f49…` on the genome-wide sample;
[details](validation/2026-10-15-subanchor-outside-repeats/README.md)). The whole-genome section was
re-run with v0.1.5; the other sections are v0.1.1.

## Whole genomes, reads to sorted and indexed BAM

ONT's GIAB 2025.01 data, R10.4.1, dorado sup v5.0.0. HG005: all reads of flow cell PAW87816
(27.9x); HG002: 60% of PAW71238 (27.2x). minimap2 `-ax lr:hq | samtools sort -@16 --write-index`
against `bonefish -o sample.bam` (built-in sort and index). Clair3 2.0.3 (sup v5.0.0 model) and
rtg vcfeval against GIAB v4.2.1; Sniffles2 and Truvari against GIAB v5.0q (HG002). bonefish v0.1.5
(15 Oct; minimap2 results from 12 Oct). Per-sample reports and raw results:
[validation/2026-10-15-giab-wgs-v2](validation/2026-10-15-giab-wgs-v2/README.md) (v0.1.1:
[validation/2026-10-12-giab-wgs](validation/2026-10-12-giab-wgs/README.md)).

| | HG005 minimap2 | HG005 bonefish | HG002 minimap2 | HG002 bonefish |
|---|---:|---:|---:|---:|
| reads to sorted, indexed BAM | 114.7 min | **34.8 min (3.3x)*** | 114.8 min | **38.4 min (3.0x)** |
| CPU time | 25.3 h | 6.9 h | 25.3 h | 6.8 h |
| peak memory (process) | 16.7 GB + sort 7.3 GB | 19.9 GB | 16.9 GB + sort 7.2 GB | 19.9 GB |
| SNP F1 | 0.9975 | 0.9975 | 0.9971 | 0.9971 |
| indel F1 | 0.9395 | 0.9395 | 0.8737 | 0.8737 |
| SV F1 vs v5.0q, strict / multi-match | – | – | 0.8055 / 0.8706 | 0.8053 / 0.8706 |
| primary alignments identical / moved | | 90.6% / 1.1% | | 95.0% / 1.4% |
| reads newly unmapped | | 3 | | 2 |

\*v0.1.1 timing: the v0.1.5 run of HG005 shared its NTFS drive with a download, which stalled the
sort's final write (56 min; its CPU time, 6.9 h, is unaffected). v0.1.1 had a higher SV F1 (0.8094 /
0.8741), from alignments through the wrong repeat copy in tandem repeats that Truvari's position-and-size
matching counts as true positives; the tandem-repeat benchmark, which realigns haplotypes, scores
v0.1.5 higher ([details](validation/2026-10-15-giab-wgs-v2/README.md#structural-variants-and-tandem-repeats)).

HG002 small variants re-scored against GIAB v5.0q (the T2T-Q100-based truth, 2.74 Gbp including
chrX/Y and difficult regions):

| region | SNP F1 minimap2 (FP) | SNP F1 bonefish (FP) | indel F1 minimap2 | indel F1 bonefish |
|---|---:|---:|---:|---:|
| all | 0.9862 (46,658) | 0.9863 (45,748) | 0.7075 | 0.7075 |
| new vs v4.2.1 | 0.8774 (42,180) | 0.8783 (41,392) | 0.4475 | 0.4475 |
| segmental duplications | 0.9333 (8,857) | 0.9335 (8,833) | 0.6506 | 0.6509 |
| not in difficult regions | 0.9971 | 0.9972 | 0.9584 | 0.9586 |

## Methylation calls

HG002 chr20 (55x, dorado sup v5.0.0 with 5mC/5hmC models), aligned by both with `-y`, then
`modkit pileup --cpg --combine-strands`: per-CpG 5mC agrees with minimap2's at r = 0.9988 (sites
>= 10x in both), mean absolute difference 0.21 percentage points, chr20 5mC level 64.42% vs 64.38%;
the few discordant sites are where the aligners placed different reads. Details:
[validation/2026-10-14-modkit-methylation](validation/2026-10-14-modkit-methylation/README.md).

## ONT's wf-human-variation pipeline

The same HG002 chr20 reads aligned by minimap2 and by bonefish, then ONT's
[wf-human-variation](https://github.com/epi2me-labs/wf-human-variation) v2.8.1 with every analysis on
(Clair3, Sniffles2, modkit, Straglr, Spectre), unchanged:

| | minimap2 BAM | bonefish BAM |
|---|---:|---:|
| SNP F1 (vs GIAB v4.2.1) | 0.9987 | 0.9987 (identical TP/FP/FN) |
| indel F1 | 0.9014 | 0.9018 |
| SV F1 (vs GIAB v5.0q) | 0.7817 | 0.7879 |
| CpG 5mC level / per-site r | 64.62% | 64.66% / 0.9988 |
| STR loci (chr20 catalogue), CNV calls | 2, 1 DUP | 2 (same alleles), same DUP |

Details: [validation/2026-10-14-wf-human-variation](validation/2026-10-14-wf-human-variation/README.md).

## Other references

Same genome-wide HG001 sample, bonefish vs minimap2: on T2T-CHM13 41 s vs 166 s (4.0x), both mapping
96.01% of reads; on GRCh38 with 261 ALT contigs and decoys (`--alt` given to both) 49 s vs 206 s (4.2x),
with the same number of primaries on ALT contigs (155 vs 154) and the same MHC mapping quality.
Details: [validation/2026-10-14-other-references](validation/2026-10-14-other-references/README.md).

## Simulated reads with known origin

Badread reads from every GRCh38 chromosome; a primary alignment is wrong if it misses the read's true
interval. bonefish against minimap2 `lr:hq` (wrong placements among reads mapped at MAPQ >= 10):
R10-like reads (65,327 scored) 0 vs 0; ultra-long (7,475) 0 vs 0; 0.5 kb reads (190,294) 19 vs 19;
93%-identity reads (31,940) 3 vs 2, where minimap2 `map-ont` has 0 and maps 257 more reads.
Details: [validation/2026-10-14-simulated-reads](validation/2026-10-14-simulated-reads/README.md).

## HAC, R9.4.1 and phasing

HG002 chr20 reads, Clair3 and Sniffles2 (minimap2 / bonefish F1): R10.4.1 sup SNP 0.9987 / 0.9987,
indel 0.9005 / 0.9013, SV 0.8011 / 0.8017; HAC SNP 0.9985 / 0.9984, indel 0.8294 / 0.8310, SV
0.7887 / 0.7917; R9.4.1 SNP 0.9980 / 0.9980, indel 0.8173 / 0.8172, SV 0.7704 / 0.7734. Phasing with
longphase against the Q100 truth: 7 switch errors on chr20 from either aligner's BAM, block N50 2.17 /
2.20 Mb. Details: [validation/2026-10-14-hac-r9-phasing](validation/2026-10-14-hac-r9-phasing/README.md).

## Medically relevant genes and tandem repeats

HG002 whole-genome calls (minimap2 / bonefish v0.1.5 F1). GIAB CMRG v1.00: SNPs 0.9650 / 0.9650,
indels 0.7514 / 0.7518, SVs 0.9082 / 0.9070. GIAB tandem-repeat benchmark v1.0.1, Tier 1, after
`truvari refine`: 0.9026 / 0.9023 (v0.1.1-v0.1.3: 0.9000, from misaligned reads in VNTRs; see
[validation/2026-10-15-tandem-repeat-detours](validation/2026-10-15-tandem-repeat-detours/README.md)).
Details: [validation/2026-10-15-giab-wgs-v2](validation/2026-10-15-giab-wgs-v2/README.md) and
[validation/2026-10-14-cmrg-tandem-repeats](validation/2026-10-14-cmrg-tandem-repeats/README.md).

## Output checks and limits

Sorted BAMs from FASTQ and from dorado BAMs with methylation tags pass Picard ValidateSamFile (no
errors with the reference), with NM/MD equal to samtools' recomputation and valid SA, RG and MM/ML
tags. Output is identical across thread counts, batch sizes, input forms (including stdin and pipes),
index types and SIMD paths. A 4.6 Gbp reference with positions past 2^32 and a 1.75 Gbp contig gives
the same alignments, and the sorted BAM gets a CSI index. Details:
[validation/2026-10-14-output-checks](validation/2026-10-14-output-checks/README.md).

## Threads and memory

Genome-wide HG001 sample: 360 s on 1 thread, 103 s on 4, 51 s on 16 (minimap2: 385 s on 4, 191 s on
16), identical output at every count. Peak memory at 16 threads is 12.3 GB with the default batch
(`-K 100M`), and 8.6 GB with `-t 4 -K 10M`. Details:
[validation/2026-10-14-scaling-memory](validation/2026-10-14-scaling-memory/README.md).

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
- Segmental duplications against GIAB v4.2.1: about 100-130 more SNP false positives per genome
  (v0.1.5; F1 −0.0005); against v5.0q bonefish has fewer false positives there.
- About 0.4% fewer read bases are aligned: low-identity read ends are left soft-clipped (primary
  alignment identity is higher, 0.979 vs 0.977).
- Clair3 takes 9–15% longer on bonefish BAMs.
- Peak memory for a human genome to sorted BAM is about 20 GB with 16 threads (measured with the
  old 500 Mbp batch; the 100 Mbp default lowers it by about 1.5 GB). The GRCh38 index is 7.3 GB, or
  8.5 GB in the raw format that loads in under a second.
- Up to v0.1.3, SV calls in tandem repeats had more false positives (+463 of 50 bp or more in GIAB's
  Tier 1 repeats; refined F1 0.9000 vs 0.9026). Fixed from v0.1.4: genome-wide F1 0.9023 (v0.1.5).
- For low-accuracy reads (simulated, ~93% identity), minimap2 `map-ont` maps about 1% more reads.
- Alignment heuristics depend slightly on absolute reference positions, as minimap2's do: shifting a
  chromosome by 1 kb changes about 3% of the reads on it (minimap2: the same).

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
