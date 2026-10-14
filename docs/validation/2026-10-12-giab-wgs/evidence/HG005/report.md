# HG005: bonefish vs minimap2, whole genome

## Input

- Reads: 6,928,131; primary query bases 86.4 Gbp (~27.9x, unmapped reads' bases not counted)
- Download of ONT's BAM from S3: 34.4 min; BAM -> FASTQ: 21.7 min

## Speed and resources

| aligner | reads -> sorted BAM | aligner wall | CPU (aligner + sort) | aligner peak RSS | sort peak RSS | system mem peak | avg CPU | avg iowait | BAM size |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| minimap2 | 114.7 min | 104.1 min | 25.29 h | 16.7 GB | 7.3 GB | 28.7 GB | 86% | 0.3% | 63.2 GB |
| bonefish | 34.8 min | 34.7 min | 7.23 h | 19.6 GB | (built in) | 18.5 GB | 82% | 0.8% | 63.3 GB |

bonefish is **3.29x** faster than minimap2 + samtools sort from reads to sorted, indexed BAM.

### Every stage (from the 10 s resource samples)

System memory includes other programs running on the machine; swap growth is relative to the swap in use when the run started (10.5 GB).

| stage | duration | avg CPU | avg iowait | system mem peak | largest process RSS | swap growth |
|---|---:|---:|---:|---:|---:|---:|
| truth | 0.7 min | 9% | 0.0% | 7.5 GB | 0.6 GB | 0.0 GB |
| download | 37.8 min | 5% | 0.3% | 7.8 GB | 0.2 GB | 0.0 GB |
| convert | 21.7 min | 46% | 0.7% | 8.0 GB | 0.0 GB | 0.0 GB |
| align_minimap2 | 112.3 min | 86% | 0.3% | 28.7 GB | 15.8 GB | 10.2 GB |
| bamstats_minimap2 | 4.7 min | 26% | 1.1% | 6.5 GB | 0.2 GB | 0.4 GB |
| clair3_minimap2 | 74.2 min | 75% | 0.2% | 20.7 GB | 6.0 GB | 3.8 GB |
| vcfeval_minimap2 | 4.0 min | 28% | 0.0% | 25.3 GB | 19.9 GB | 0.6 GB |
| sniffles_minimap2 | 1.7 min | 67% | 8.9% | 12.1 GB | 1.3 GB | 0.9 GB |
| align_bonefish | 34.2 min | 82% | 0.8% | 18.5 GB | 19.4 GB | 3.1 GB |
| bamstats_bonefish | 4.5 min | 27% | 1.4% | 5.5 GB | 0.3 GB | 1.0 GB |
| clair3_bonefish | 80.7 min | 74% | 1.2% | 28.2 GB | 4.8 GB | 15.3 GB |
| vcfeval_bonefish | 4.0 min | 28% | 0.0% | 26.6 GB | 20.7 GB | 1.7 GB |
| sniffles_bonefish | 1.7 min | 66% | 9.6% | 13.9 GB | 1.4 GB | 1.8 GB |
| compare | 0.3 min | 6% | 0.5% | 6.4 GB | 1.5 GB | 1.8 GB |

## BAM statistics

| statistic | minimap2 | bonefish |
|---|---:|---:|
| reads | 6928131 | 6928131 |
| primary_mapped | 6431077 | 6431077 |
| unmapped | 497054 | 497054 |
| pct_reads_mapped | 92.826 | 92.826 |
| supplementary | 636262 | 663472 |
| secondary | 1781816 | 1649949 |
| pct_primary_bases_aligned | 96.581 | 96.210 |
| primary_identity | 0.97679 | 0.97885 |
| pct_primary_mapq_ge_20 | 94.009 | 94.136 |
| mapq_0-0 | 39466 | 39436 |
| mapq_1-9 | 288138 | 279693 |
| mapq_10-19 | 57692 | 58000 |
| mapq_20-29 | 45622 | 45474 |
| mapq_30-59 | 131873 | 133214 |
| mapq_60-60 | 5868286 | 5875260 |

### Primary alignments, bonefish vs minimap2

| category | reads |
|---|---:|
| compared | 6928131 |
| identical | 6275403 |
| same_locus | 76360 |
| moved | 79311 |
| both_unmapped | 497051 |
| newly_unmapped | 3 |
| newly_mapped | 3 |
| mapq_lower | 28064 |
| mapq_higher | 48860 |
| mapq_lost_20 | 7286 |
| mapq_gained_20 | 15453 |
| only_in_test | 0 |
| only_in_ref | 0 |

## Clair3 (sup v5.0.0 model, longphase intermediate phasing)

| BAM from | wall | CPU | system mem peak | calls | PASS SNPs | PASS indels |
|---|---:|---:|---:|---:|---:|---:|
| minimap2 | 74.7 min | 14.39 h | 20.7 GB | 5409848 | 4361151 | 898785 |
| bonefish | 81.3 min | 15.42 h | 28.2 GB | 5482463 | 4429712 | 895350 |

### Small-variant accuracy vs GIAB v4.2.1 (rtg vcfeval, PASS calls)

| region / type | minimap2 F1 | TP/FP/FN | bonefish F1 | TP/FP/FN | ΔF1 |
|---|---:|---:|---:|---:|---:|
| all indels | 0.9395 | 392231/26098/24434 | 0.9396 | 392155/25948/24510 | +0.0001 |
| all snps | 0.9975 | 3268693/8774/7778 | 0.9974 | 3268656/8915/7815 | -0.0001 |
| AllTandemRepeatsandHomopolymers_slop5 indels | 0.8965 | 210068/25217/23293 | 0.8967 | 210010/25064/23351 | +0.0002 |
| AllTandemRepeatsandHomopolymers_slop5 snps | 0.9901 | 148018/1909/1053 | 0.9901 | 148015/1909/1056 | +0.0000 |
| alldifficultregions indels | 0.9064 | 238284/25559/23656 | 0.9065 | 238224/25415/23716 | +0.0001 |
| alldifficultregions snps | 0.9910 | 572721/6733/3651 | 0.9909 | 572708/6876/3664 | -0.0001 |
| lowmappabilityall indels | 0.9728 | 9500/242/290 | 0.9722 | 9499/253/291 | -0.0006 |
| lowmappabilityall snps | 0.9855 | 190880/3540/2091 | 0.9850 | 190858/3693/2113 | -0.0005 |
| notinalldifficultregions indels | 0.9957 | 154067/536/788 | 0.9957 | 154051/530/804 | +0.0000 |
| notinalldifficultregions snps | 0.9989 | 2695972/2041/4127 | 0.9989 | 2695948/2039/4151 | +0.0000 |
| segdups indels | 0.9660 | 10039/392/315 | 0.9653 | 10038/406/316 | -0.0007 |
| segdups snps | 0.9735 | 122070/4594/2056 | 0.9728 | 122048/4750/2078 | -0.0007 |

## Structural variants (Sniffles2)

| BAM from | Sniffles wall | F1 vs GIAB v5.0q, strict (TP/FP/FN) | F1, one call may match several truth SVs |
|---|---:|---:|---:|
| minimap2 | 1.7 min | - | - |
| bonefish | 1.7 min | - | - |

SV call agreement (Truvari, minimap2 calls as base): minimap2_svs_found_by_bonefish 24063; bonefish_only 1802; minimap2_only 1272; agreement_f1 0.9400

