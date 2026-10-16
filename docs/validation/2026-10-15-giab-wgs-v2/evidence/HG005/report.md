# HG005_v2: bonefish vs minimap2, whole genome

## Input

- Reads: 6,928,131; primary query bases 86.4 Gbp (~27.9x, unmapped reads' bases not counted)

## Speed and resources

| aligner | reads -> sorted BAM | aligner wall | CPU (aligner + sort) | aligner peak RSS | sort peak RSS | system mem peak | avg CPU | avg iowait | BAM size |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| minimap2 | 114.7 min | 104.1 min | 25.29 h | 16.7 GB | 7.3 GB | - | - | - | 63.2 GB |
| bonefish | 56.4 min | 56.3 min | 6.95 h | 19.9 GB | (built in) | 17.7 GB | 57% | 4.8% | 63.3 GB |

bonefish is **2.03x** faster than minimap2 + samtools sort from reads to sorted, indexed BAM.

### Every stage (from the 10 s resource samples)

System memory includes other programs running on the machine; swap growth is relative to the swap in use when the run started (5.0 GB).

| stage | duration | avg CPU | avg iowait | system mem peak | largest process RSS | swap growth |
|---|---:|---:|---:|---:|---:|---:|
| align_bonefish | 51.7 min | 57% | 4.8% | 17.7 GB | 19.8 GB | 4.0 GB |
| bamstats_bonefish | 5.2 min | 30% | 1.7% | 4.6 GB | 0.3 GB | 0.7 GB |
| clair3_bonefish | 96.5 min | 66% | 11.7% | 30.0 GB | 4.6 GB | 14.8 GB |
| vcfeval_bonefish | 3.8 min | 28% | 0.0% | 24.5 GB | 20.5 GB | 0.9 GB |
| sniffles_bonefish | 1.7 min | 66% | 9.6% | 12.7 GB | 1.8 GB | 1.0 GB |
| compare | 0.3 min | 7% | 0.5% | 5.6 GB | 1.2 GB | 1.0 GB |

## BAM statistics

| statistic | minimap2 | bonefish |
|---|---:|---:|
| reads | 6928131 | 6928131 |
| primary_mapped | 6431077 | 6431077 |
| unmapped | 497054 | 497054 |
| pct_reads_mapped | 92.826 | 92.826 |
| supplementary | 636262 | 662671 |
| secondary | 1781816 | 1649598 |
| pct_primary_bases_aligned | 96.581 | 96.212 |
| primary_identity | 0.97679 | 0.97896 |
| pct_primary_mapq_ge_20 | 94.009 | 94.138 |
| mapq_0-0 | 39466 | 39422 |
| mapq_1-9 | 288138 | 279689 |
| mapq_10-19 | 57692 | 57894 |
| mapq_20-29 | 45622 | 45521 |
| mapq_30-59 | 131873 | 133155 |
| mapq_60-60 | 5868286 | 5875396 |

### Primary alignments, bonefish vs minimap2

| category | reads |
|---|---:|
| compared | 6928131 |
| identical | 6278002 |
| same_locus | 75282 |
| moved | 77790 |
| both_unmapped | 497051 |
| newly_unmapped | 3 |
| newly_mapped | 3 |
| mapq_lower | 23482 |
| mapq_higher | 45115 |
| mapq_lost_20 | 6306 |
| mapq_gained_20 | 14597 |
| only_in_test | 0 |
| only_in_ref | 0 |

## Clair3 (sup v5.0.0 model, longphase intermediate phasing)

| BAM from | wall | CPU | system mem peak | calls | PASS SNPs | PASS indels |
|---|---:|---:|---:|---:|---:|---:|
| minimap2 | 74.7 min | 14.39 h | - | 5409848 | 4361151 | 898785 |
| bonefish | 98.2 min | 15.72 h | 30.0 GB | 5499001 | 4441946 | 899694 |

### Small-variant accuracy vs GIAB v4.2.1 (rtg vcfeval, PASS calls)

| region / type | minimap2 F1 | TP/FP/FN | bonefish F1 | TP/FP/FN | ΔF1 |
|---|---:|---:|---:|---:|---:|
| all indels | 0.9395 | 392231/26098/24434 | 0.9395 | 392158/25958/24507 | +0.0000 |
| all snps | 0.9975 | 3268693/8774/7778 | 0.9975 | 3268660/8887/7811 | +0.0000 |
| AllTandemRepeatsandHomopolymers_slop5 indels | 0.8965 | 210068/25217/23293 | 0.8966 | 210014/25080/23347 | +0.0001 |
| AllTandemRepeatsandHomopolymers_slop5 snps | 0.9901 | 148018/1909/1053 | 0.9901 | 148019/1904/1052 | +0.0000 |
| alldifficultregions indels | 0.9064 | 238284/25559/23656 | 0.9065 | 238226/25425/23714 | +0.0001 |
| alldifficultregions snps | 0.9910 | 572721/6733/3651 | 0.9909 | 572714/6852/3658 | -0.0001 |
| lowmappabilityall indels | 0.9728 | 9500/242/290 | 0.9724 | 9497/246/293 | -0.0004 |
| lowmappabilityall snps | 0.9855 | 190880/3540/2091 | 0.9851 | 190859/3673/2112 | -0.0004 |
| notinalldifficultregions indels | 0.9957 | 154067/536/788 | 0.9957 | 154052/530/803 | +0.0000 |
| notinalldifficultregions snps | 0.9989 | 2695972/2041/4127 | 0.9989 | 2695946/2035/4153 | +0.0000 |
| segdups indels | 0.9660 | 10039/392/315 | 0.9656 | 10036/397/318 | -0.0004 |
| segdups snps | 0.9735 | 122070/4594/2056 | 0.9729 | 122045/4725/2081 | -0.0006 |

## Structural variants (Sniffles2)

| BAM from | Sniffles wall | F1 vs GIAB v5.0q, strict (TP/FP/FN) | F1, one call may match several truth SVs |
|---|---:|---:|---:|
| minimap2 | 1.7 min | - | - |
| bonefish | 1.7 min | - | - |

SV call agreement (Truvari, minimap2 calls as base): minimap2_svs_found_by_bonefish 24778; bonefish_only 496; minimap2_only 557; agreement_f1 0.9792

