# HG002_v2: bonefish vs minimap2, whole genome

## Input

- Reads: 4,751,476; primary query bases 84.3 Gbp (~27.2x, unmapped reads' bases not counted)

## Speed and resources

| aligner | reads -> sorted BAM | aligner wall | CPU (aligner + sort) | aligner peak RSS | sort peak RSS | system mem peak | avg CPU | avg iowait | BAM size |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| minimap2 | 114.8 min | 104.8 min | 25.27 h | 16.9 GB | 7.2 GB | - | - | - | 59.8 GB |
| bonefish | 38.4 min | 38.3 min | 6.76 h | 19.9 GB | (built in) | 16.2 GB | 68% | 4.2% | 59.8 GB |

bonefish is **2.99x** faster than minimap2 + samtools sort from reads to sorted, indexed BAM.

### Every stage (from the 10 s resource samples)

System memory includes other programs running on the machine; swap growth is relative to the swap in use when the run started (6.0 GB).

| stage | duration | avg CPU | avg iowait | system mem peak | largest process RSS | swap growth |
|---|---:|---:|---:|---:|---:|---:|
| align_bonefish | 37.5 min | 68% | 4.2% | 16.2 GB | 19.7 GB | 1.2 GB |
| bamstats_bonefish | 4.8 min | 24% | 1.5% | 4.1 GB | 0.3 GB | 0.2 GB |
| clair3_bonefish | 91.5 min | 69% | 5.4% | 28.8 GB | 4.3 GB | 11.9 GB |
| vcfeval_bonefish | 4.0 min | 27% | 0.0% | 24.5 GB | 20.3 GB | 0.2 GB |
| sniffles_bonefish | 1.7 min | 66% | 7.6% | 12.4 GB | 1.0 GB | 0.2 GB |
| truvari_bonefish | 4.0 min | 6% | 0.0% | 4.3 GB | 0.0 GB | 0.2 GB |
| compare | 0.2 min | 6% | 0.0% | 5.0 GB | 0.8 GB | 0.2 GB |

## BAM statistics

| statistic | minimap2 | bonefish |
|---|---:|---:|
| reads | 4751476 | 4751476 |
| primary_mapped | 4647266 | 4647264 |
| unmapped | 104210 | 104212 |
| pct_reads_mapped | 97.807 | 97.807 |
| supplementary | 591184 | 605858 |
| secondary | 1364104 | 1250799 |
| pct_primary_bases_aligned | 95.701 | 95.304 |
| primary_identity | 0.97544 | 0.97742 |
| pct_primary_mapq_ge_20 | 94.373 | 94.558 |
| mapq_0-0 | 16074 | 16057 |
| mapq_1-9 | 213867 | 205348 |
| mapq_10-19 | 31540 | 31522 |
| mapq_20-29 | 24499 | 24163 |
| mapq_30-59 | 66358 | 67346 |
| mapq_60-60 | 4294928 | 4302828 |

### Primary alignments, bonefish vs minimap2

| category | reads |
|---|---:|
| compared | 4751476 |
| identical | 4512729 |
| same_locus | 68479 |
| moved | 66056 |
| both_unmapped | 104210 |
| newly_unmapped | 2 |
| newly_mapped | 0 |
| mapq_lower | 20371 |
| mapq_higher | 39833 |
| mapq_lost_20 | 5420 |
| mapq_gained_20 | 13972 |
| only_in_test | 0 |
| only_in_ref | 0 |

## Clair3 (sup v5.0.0 model, longphase intermediate phasing)

| BAM from | wall | CPU | system mem peak | calls | PASS SNPs | PASS indels |
|---|---:|---:|---:|---:|---:|---:|
| minimap2 | 80.4 min | 15.43 h | - | 5466877 | 4399682 | 912324 |
| bonefish | 92.8 min | 16.53 h | 28.8 GB | 5585712 | 4505957 | 914387 |

### Small-variant accuracy vs GIAB v4.2.1 (rtg vcfeval, PASS calls)

| region / type | minimap2 F1 | TP/FP/FN | bonefish F1 | TP/FP/FN | ΔF1 |
|---|---:|---:|---:|---:|---:|
| all indels | 0.8737 | 457142/63973/68255 | 0.8737 | 457105/63860/68292 | +0.0000 |
| all snps | 0.9971 | 3358537/11752/7837 | 0.9971 | 3358542/11785/7832 | +0.0000 |
| AllTandemRepeatsandHomopolymers_slop5 indels | 0.8069 | 271816/62892/67238 | 0.8070 | 271792/62789/67262 | +0.0001 |
| AllTandemRepeatsandHomopolymers_slop5 snps | 0.9850 | 184833/3506/2107 | 0.9851 | 184831/3478/2109 | +0.0001 |
| alldifficultregions indels | 0.8208 | 300191/63450/67601 | 0.8209 | 300162/63340/67630 | +0.0001 |
| alldifficultregions snps | 0.9875 | 609467/10632/4812 | 0.9874 | 609457/10697/4822 | -0.0001 |
| lowmappabilityall indels | 0.9633 | 10117/481/290 | 0.9642 | 10116/460/291 | +0.0009 |
| lowmappabilityall snps | 0.9812 | 190389/5062/2255 | 0.9815 | 190382/4899/2262 | +0.0003 |
| notinalldifficultregions indels | 0.9962 | 157125/521/664 | 0.9962 | 157119/519/670 | +0.0000 |
| notinalldifficultregions snps | 0.9992 | 2749070/1120/3025 | 0.9993 | 2749085/1088/3010 | +0.0001 |
| segdups indels | 0.9514 | 10475/718/353 | 0.9519 | 10472/703/356 | +0.0005 |
| segdups snps | 0.9618 | 118741/7067/2374 | 0.9613 | 118731/7171/2384 | -0.0005 |

## Structural variants (Sniffles2)

| BAM from | Sniffles wall | F1 vs GIAB v5.0q, strict (TP/FP/FN) | F1, one call may match several truth SVs |
|---|---:|---:|---:|
| minimap2 | 1.6 min | 0.8055 (20127/1727/7996) | 0.8706 (23085/1601/5038) |
| bonefish | 1.6 min | 0.8053 (20096/1693/8027) | 0.8706 (23055/1569/5068) |

SV call agreement (Truvari, minimap2 calls as base): minimap2_svs_found_by_bonefish 24973; bonefish_only 550; minimap2_only 590; agreement_f1 0.9777

