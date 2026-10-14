# ONT's wf-human-variation on bonefish vs minimap2 BAMs (14 Oct 2026)

ONT's own human variation pipeline, [epi2me-labs/wf-human-variation](https://github.com/epi2me-labs/wf-human-variation)
v2.8.1, run unchanged on two BAMs of the same reads: one aligned by minimap2, one by bonefish. Every
analysis the pipeline offers was switched on: small variants (Clair3), structural variants (Sniffles2),
methylation (modkit), short tandem repeats (Straglr) and copy-number variants (Spectre).

## Data and setup

- **Reads:** every primary read on chr20 of ONT's HG002 flow cell PAW70337 (GIAB 2025.01, R10.4.1,
  dorado sup v5.0.0 with 5mC/5hmC and 6mA models), 188,140 reads, 3.52 Gbp (~55x of chr20), with
  `MM`/`ML` tags. These are the BAMs of the [methylation comparison](../2026-10-14-modkit-methylation/README.md):
  `minimap2 -ax lr:hq -y | samtools sort` and `bonefish -y -o`, both against the whole GRCh38 no-alt genome.
- **Pipeline:** `nextflow run epi2me-labs/wf-human-variation -r v2.8.1 --bam X.bam --ref GRCh38_no_alt.fa
  --bed chr20.bed --snp --sv --mod --str --cnv --sex XY --override_basecaller_cfg
  dna_r10.4.1_e8.2_400bps_sup@v5.0.0`, Docker, 16 CPUs and 26 GB for the local executor
  (`evidence/local.config`). Nextflow 26.10 with `NXF_SYNTAX_PARSER=v1` (v2.8.1 predates the strict
  parser). The basecall model is passed explicitly because the FASTQ-derived BAMs have no `@RG` lines.
- **Scoring** (`evidence/wfhv_score.sh`): PASS SNPs and indels with rtg vcfeval against GIAB v4.2.1
  (chr20 benchmark regions); SVs with Truvari against GIAB v5.0q (chr20); methylation and STRs compared
  between the two runs; CNV calls listed.

## Results

| analysis | from minimap2's BAM | from bonefish's BAM |
|---|---:|---:|
| SNP F1 (TP / FP / FN) | 0.9987 (71,267 / 90 / 95) | 0.9987 (71,267 / 90 / 95) |
| indel F1 (TP / FP / FN) | 0.9014 (9,914 / 828 / 1,342) | 0.9018 (9,921 / 826 / 1,335) |
| SV F1 vs v5.0q (TP / FP / FN) | 0.7817 (548 / 52 / 254) | 0.7879 (559 / 58 / 243) |
| CpG 5mC level (chr20) | 64.62% | 64.66% |
| per-CpG 5mC agreement (>= 10x in both) | r = 0.9988 (751,483 sites) | |
| STR genotypes (Straglr, chr20 catalogue loci) | 2 | 2, identical alleles |
| CNV calls (Spectre) | DUP chr20:34,421,000-34,531,000, CN 3 | DUP chr20:34,420,000-34,531,000, CN 3 |
| pipeline wall time | 18 min 16 s | 20 min 31 s |

- SNP calls are identical in count, and indel and SV accuracy are slightly higher from bonefish's BAM.
- The methylation agreement matches the direct modkit comparison: the few discordant CpGs are where the
  aligners placed different reads.
- The CNV call is the same duplication; its start differs by one 1 kb bin.
- The pipeline ran 2 min longer on bonefish's BAM, mostly in Clair3 (30.0 vs 27.7 task-minutes in the
  Nextflow trace), which also takes longer on bonefish BAMs in the whole-genome runs
  ([GIAB whole genomes](../2026-10-12-giab-wgs/README.md)).
  bonefish aligned the reads to a sorted BAM in under 2 minutes against 8 for minimap2 + samtools.

## Limits

One chromosome, one sample. The chr20 STR catalogue holds only two loci, so the STR comparison shows the
pipeline runs end to end on bonefish BAMs rather than measuring STR accuracy. The pipeline's own alignment
step was not used: it aligns unaligned BAMs with minimap2, so the aligned BAMs were given as input.

## Files

- `evidence/score/`: the comparison tables (`small_variants.tsv`, `sv.tsv`, `mods.tsv`, `str.tsv`,
  `cnv.txt`) and the rtg and Truvari summaries.
- `evidence/<aligner>/`: the pipeline's SV, CNV and STR calls, read statistics, coverage summary,
  `params.json`, the Nextflow task trace and wall time.
- `evidence/wfhv_compare.sh`, `wfhv_score.sh`, `local.config`, `chr20.bed`, `run.log`.
- `evidence/OUTPUTS_NOT_IN_GIT.sha256`: SHA-256 of the SNP VCFs, bedMethyl, bigWig and CRAM outputs.
- `SHA256SUMS` covers `evidence/`. Local paths are replaced by `${REPO}`, `${DATA}` and `${HOME}`.
