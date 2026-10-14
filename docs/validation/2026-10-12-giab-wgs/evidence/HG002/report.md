# HG002: bonefish vs minimap2, whole genome

## Input

- Reads: 4,751,476; primary query bases 84.3 Gbp (~27.2x, unmapped reads' bases not counted)
- Download of ONT's BAM from S3: 54.4 min; BAM -> FASTQ: 22.0 min

## Speed and resources

| aligner | reads -> sorted BAM | aligner wall | CPU (aligner + sort) | aligner peak RSS | sort peak RSS | system mem peak | avg CPU | avg iowait | BAM size |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| minimap2 | 114.8 min | 104.8 min | 25.27 h | 16.9 GB | 7.2 GB | 29.3 GB | 85% | 0.2% | 59.8 GB |
| bonefish | 37.5 min | 37.4 min | 7.17 h | 19.7 GB | (built in) | 17.9 GB | 80% | 1.5% | 59.8 GB |

bonefish is **3.06x** faster than minimap2 + samtools sort from reads to sorted, indexed BAM.

### Every stage (from the 10 s resource samples)

System memory includes other programs running on the machine; swap growth is relative to the swap in use when the run started (12.3 GB).

| stage | duration | avg CPU | avg iowait | system mem peak | largest process RSS | swap growth |
|---|---:|---:|---:|---:|---:|---:|
| truth | 0.7 min | 8% | 0.0% | 6.0 GB | 0.6 GB | 0.0 GB |
| download | 54.2 min | 3% | 0.1% | 7.2 GB | 0.2 GB | 0.0 GB |
| convert | 21.8 min | 49% | 0.7% | 7.0 GB | 0.0 GB | 0.0 GB |
| align_minimap2 | 112.7 min | 85% | 0.2% | 29.3 GB | 16.5 GB | 4.4 GB |
| bamstats_minimap2 | 4.8 min | 25% | 0.9% | 6.0 GB | 0.3 GB | 0.0 GB |
| clair3_minimap2 | 79.5 min | 76% | 0.2% | 26.3 GB | 5.5 GB | 0.0 GB |
| vcfeval_minimap2 | 4.0 min | 28% | 0.1% | 26.9 GB | 20.7 GB | 0.0 GB |
| sniffles_minimap2 | 1.7 min | 70% | 7.2% | 14.0 GB | 1.0 GB | 0.0 GB |
| truvari_minimap2 | 4.0 min | 7% | 0.0% | 6.5 GB | 0.0 GB | 0.0 GB |
| align_bonefish | 36.8 min | 80% | 1.5% | 17.9 GB | 19.7 GB | 2.4 GB |
| bamstats_bonefish | 5.2 min | 31% | 1.4% | 7.1 GB | 0.3 GB | 0.0 GB |
| clair3_bonefish | 89.0 min | 76% | 0.6% | 28.7 GB | 4.0 GB | 12.1 GB |
| vcfeval_bonefish | 4.2 min | 29% | 0.1% | 27.1 GB | 20.7 GB | 0.4 GB |
| sniffles_bonefish | 1.7 min | 68% | 8.9% | 12.8 GB | 1.1 GB | 1.8 GB |
| truvari_bonefish | 4.2 min | 11% | 0.6% | 7.3 GB | 0.0 GB | 0.6 GB |
| compare | 0.3 min | 11% | 0.0% | 7.6 GB | 0.5 GB | 0.2 GB |

## BAM statistics

| statistic | minimap2 | bonefish |
|---|---:|---:|
| reads | 4751476 | 4751476 |
| primary_mapped | 4647266 | 4647264 |
| unmapped | 104210 | 104212 |
| pct_reads_mapped | 97.807 | 97.807 |
| supplementary | 591184 | 606684 |
| secondary | 1364104 | 1251138 |
| pct_primary_bases_aligned | 95.701 | 95.302 |
| primary_identity | 0.97544 | 0.97732 |
| pct_primary_mapq_ge_20 | 94.373 | 94.558 |
| mapq_0-0 | 16074 | 16076 |
| mapq_1-9 | 213867 | 205351 |
| mapq_10-19 | 31540 | 31489 |
| mapq_20-29 | 24499 | 24247 |
| mapq_30-59 | 66358 | 67477 |
| mapq_60-60 | 4294928 | 4302624 |

### Primary alignments, bonefish vs minimap2

| category | reads |
|---|---:|
| compared | 4751476 |
| identical | 4510433 |
| same_locus | 69738 |
| moved | 67093 |
| both_unmapped | 104210 |
| newly_unmapped | 2 |
| newly_mapped | 0 |
| mapq_lower | 24182 |
| mapq_higher | 42697 |
| mapq_lost_20 | 6135 |
| mapq_gained_20 | 14698 |
| only_in_test | 0 |
| only_in_ref | 0 |

## Clair3 (sup v5.0.0 model, longphase intermediate phasing)

| BAM from | wall | CPU | system mem peak | calls | PASS SNPs | PASS indels |
|---|---:|---:|---:|---:|---:|---:|
| minimap2 | 80.4 min | 15.43 h | 26.3 GB | 5466877 | 4399682 | 912324 |
| bonefish | 89.7 min | 17.03 h | 28.7 GB | 5567862 | 4492319 | 910121 |

### Small-variant accuracy vs GIAB v4.2.1 (rtg vcfeval, PASS calls)

| region / type | minimap2 F1 | TP/FP/FN | bonefish F1 | TP/FP/FN | ΔF1 |
|---|---:|---:|---:|---:|---:|
| all indels | 0.8737 | 457142/63973/68255 | 0.8737 | 457092/63879/68305 | +0.0000 |
| all snps | 0.9971 | 3358537/11752/7837 | 0.9971 | 3358505/11851/7869 | +0.0000 |
| AllTandemRepeatsandHomopolymers_slop5 indels | 0.8069 | 271816/62892/67238 | 0.8069 | 271782/62808/67272 | +0.0000 |
| AllTandemRepeatsandHomopolymers_slop5 snps | 0.9850 | 184833/3506/2107 | 0.9850 | 184811/3484/2129 | +0.0000 |
| alldifficultregions indels | 0.8208 | 300191/63450/67601 | 0.8209 | 300150/63357/67642 | +0.0001 |
| alldifficultregions snps | 0.9875 | 609467/10632/4812 | 0.9873 | 609417/10759/4862 | -0.0002 |
| lowmappabilityall indels | 0.9633 | 10117/481/290 | 0.9639 | 10114/464/293 | +0.0006 |
| lowmappabilityall snps | 0.9812 | 190389/5062/2255 | 0.9813 | 190365/4978/2279 | +0.0001 |
| notinalldifficultregions indels | 0.9962 | 157125/521/664 | 0.9962 | 157118/520/671 | +0.0000 |
| notinalldifficultregions snps | 0.9992 | 2749070/1120/3025 | 0.9993 | 2749088/1092/3007 | +0.0001 |
| segdups indels | 0.9514 | 10475/718/353 | 0.9518 | 10470/703/358 | +0.0004 |
| segdups snps | 0.9618 | 118741/7067/2374 | 0.9610 | 118718/7227/2397 | -0.0008 |

## Structural variants (Sniffles2)

| BAM from | Sniffles wall | F1 vs GIAB v5.0q, strict (TP/FP/FN) | F1, one call may match several truth SVs |
|---|---:|---:|---:|
| minimap2 | 1.6 min | 0.8055 (20127/1727/7996) | 0.8706 (23085/1601/5038) |
| bonefish | 1.7 min | 0.8094 (20389/1867/7734) | 0.8741 (23363/1729/4760) |

SV call agreement (Truvari, minimap2 calls as base): minimap2_svs_found_by_bonefish 24202; bonefish_only 1869; minimap2_only 1361; agreement_f1 0.9374

