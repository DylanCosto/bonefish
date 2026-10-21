# Using bonefish in pipelines

bonefish takes the same inputs as minimap2 and writes standard SAM or BAM, so it can replace
minimap2 in a pipeline in two ways:

1. **Swap the command.** Replace the aligner call itself.
2. **Align first.** Align with bonefish, then hand the sorted BAM to a pipeline that accepts
   aligned input.

## Replacing a minimap2 command

| minimap2 | bonefish |
|---|---|
| `minimap2 -ax lr:hq -y -t 16 ref.fa reads.fq \| samtools sort -o out.bam && samtools index out.bam` | `bonefish -t 16 -y -o out.bam ref.fa reads.fq` |
| `samtools fastq -T '*' reads.bam \| minimap2 -ax lr:hq -y ...` (dorado BAM) | `bonefish -t 16 -y -o out.bam ref.fa reads.bam` |
| `minimap2 -ax splice -uf -k14 ref.fa drna.fq` | `bonefish -x rna -a ref.fa drna.fq` |
| `minimap2 -ax splice ref.fa cdna.fq` | `bonefish -x cdna -a ref.fa cdna.fq` |
| `minimap2 -ax splice --junc-bed anno.bed ...` | `bonefish -x cdna -a --junc-bed anno.bed ...` |

Notes:

- **`-o out.bam`** writes the sorted BAM and its `.bai` index in one step: the same records and
  order as `samtools sort --write-index`.
- **`-y`** keeps each read's tags (methylation `MM`/`ML`, `RG`, ...) and the input's `@RG`
  header lines. Those headers carry dorado's basecalling model, which tools such as Clair3 read.
- **`--junc-bed`** takes BED12 (one transcript per line) or BED6 (one intron per line). minimap2's
  `paftools.js gff2bed` converts a GTF or GFF3 file.
- **Indexes:** build one once per preset for repeated runs (see [usage](usage.md#indexes)).

## ONT's wf-human-variation

[wf-human-variation](https://github.com/epi2me-labs/wf-human-variation) accepts aligned BAMs.
Align with bonefish, then run the workflow on the result; it skips its own alignment step:

```bash
bonefish -t 16 -y -o sample.bam GRCh38.fa reads.bam
nextflow run epi2me-labs/wf-human-variation --bam sample.bam --ref GRCh38.fa \
    --snp --sv --mod --str --cnv --sex XY
```

- **Same reference:** use the same reference file for bonefish and the workflow.
- **Basecalling model:** with `-y` and a dorado BAM as input, the `@RG` lines carry the model. From
  FASTQ input, pass it explicitly with `--override_basecaller_cfg`, e.g.
  `dna_r10.4.1_e8.2_400bps_sup@v5.0.0`.
- **Validation:** on HG002 chr20, every analysis of the workflow (small variants, SVs,
  methylation, STRs, CNVs) gave the same or slightly better results from bonefish's BAM than from
  minimap2's ([report](validation/2026-10-14-wf-human-variation/README.md)).

## Variant callers and methylation

bonefish's sorted, indexed BAMs have been tested with these tools:

- **Small variants:** Clair3.
- **Structural variants:** Sniffles2.
- **Phasing:** longphase, WhatsHap.
- **Methylation:** modkit, which needs `-y` so the `MM`/`ML` tags are kept.

The BAMs are standard (they pass Picard ValidateSamFile), so other tools that read sorted BAMs
should work too.

Accuracy against minimap2's BAMs is in [benchmarks.md](benchmarks.md).

## RNA and cDNA

Spliced BAMs from `-x rna` or `-x cdna` are the input that genome-based transcript tools take, for example:

- **Isoform discovery and quantification:** IsoQuant, bambu, StringTie2 (`-L`), FLAIR.
- **Direct-RNA modifications:** modkit with `-y`; dorado's m6A and other calls are carried as
  `MM`/`ML` tags.

**Transcript strand:** the `ts:A` tag gives the transcript strand of each spliced alignment, as
minimap2 does.

**Not yet validated:** the presets are checked against minimap2, GENCODE junctions and simulated
reads ([report](validation/2026-10-21-rna-cdna/README.md)), but not yet against these downstream
tools.

## Containers

```bash
docker build -t bonefish https://github.com/DylanCosto/bonefish.git
docker run --rm -v "$PWD":/data bonefish bonefish -t 16 -y -o /data/out.bam /data/ref.fa /data/reads.bam
```

The image includes `procps`, which Nextflow needs for task monitoring.
