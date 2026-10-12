# bonefish

A fast aligner for Oxford Nanopore DNA reads, written in Rust. bonefish takes ONT reads (FASTQ,
or the unaligned BAM that dorado writes) and produces SAM, PAF, or a coordinate-sorted, indexed BAM
in one step. It is tuned and validated for variant calling: against GIAB truth sets, its SNP and
indel accuracy match minimap2's and its structural-variant recall is at or above minimap2's, while
aligning 2.6-8x faster.

## Quick start

```bash
cargo build --release --features jemalloc
./target/release/bonefish -d ref.raw --idx-raw ref.fa          # optional: build a fast-loading index
./target/release/bonefish -t 16 -y --sort -o sample.bam ref.raw reads.bam   # sorted BAM + .bai
./target/release/bonefish -t 16 -a ref.fa reads.fq.gz > sample.sam
```

- There is one preset, `dna` (ONT genomic DNA, R10.4.1 or R9.4.1), applied by default. The
  minimap2 names `lr:hq`, `lr:hq-fast`, `map-ont` and `lr` are accepted for it; other minimap2
  presets (PacBio, short reads, RNA/splice, assembly) are not supported.
- `-y` copies the base-modification tags (`MM`/`ML`/`MN`) from a dorado BAM (or FASTQ comments)
  into the output.
- `--sort -o FILE.bam` writes a coordinate-sorted BAM and `FILE.bam.bai` directly (same records and
  order as `bonefish -a | samtools sort --write-index`); `--sort-mem` (default 4G) bounds memory,
  with spills to `--sort-tmp`.
- Most minimap2 command-line options (`-t`, `-a`, `-c`, `--cs`, `--MD`, `--eqx`, `-N`, `-p`, ...)
  work as in minimap2.

## Speed and accuracy (16 threads)

| dataset | minimap2 `lr:hq` | bonefish | speedup |
|---|---:|---:|---:|
| HG002 R10.4.1 HAC, chr20 reads (3.5 Gbp) | 417 s | 72 s | 5.8x |
| HG002 R9.4.1, chr20 reads (4.6 Gbp) | 412 s | 103 s | 4.0x |
| HG002 R10.4.1 ultra-long (3.4 Gbp) | 542 s | 65 s | 8.3x |
| HG001 R10.4.1 sup, random genome-wide sample (2.7 Gbp) | 202 s | 45 s | 4.5x |
| *Salmonella enterica* R10.4.1 (1.2 Gbp) | 25 s | 10 s | 2.6x |

HG001 37x R10.4.1, reads to sorted indexed BAM: 40 minutes.

Variant calling (Clair3 + rtg vcfeval vs GIAB v4.2.1; Sniffles2 + Truvari vs GIAB v5.0q):

| dataset | aligner | SNP F1 | indel F1 | SV F1 |
|---|---|---:|---:|---:|
| HG002 HAC chr20 | minimap2 | 0.9985 | 0.8294 | 0.7887 |
| | bonefish | 0.9984 | 0.8310 | 0.7917 |
| HG002 R9.4.1 chr20 | minimap2 | 0.9980 | 0.8173 | 0.7704 |
| | bonefish | 0.9980 | 0.8172 | 0.7734 |

The measurements, methods and every optimisation tried (adopted or rejected) are written up in the
development repository's README.

## How it differs from minimap2

Same seeding and scoring as minimap2's `lr:hq`, with cheaper alignment and chaining: narrower,
adaptive DP bands with retries where an indel needs a wider band, a probe that stops extending
into low-identity read ends, ungapped-first and sub-anchored gap filling, a 5 kb reference-gap limit
in the main chaining pass (with a long-join rescue that skips repeat-dense reads), segment-tree
RMQ chaining and diagonal-bucketed exact chaining, plus multithreaded I/O and a built-in sorter.
Output is therefore not byte-identical to minimap2. `bonefish --help` lists the individual
options.

## Library

`bonefish-core` (imported as `bonefish`) exposes the aligner as a Rust library
(`bonefish::Aligner`, `Preset::Dna`); see `bonefish-core/examples/`.

## Credits and license

bonefish is a fork of [rammap](https://github.com/jwanglab/rammap), a pure-Rust reimplementation
of [minimap2](https://github.com/lh3/minimap2) by Heng Li. The index format, SIMD DP kernels and
most of the mapping pipeline come from rammap; the alignment algorithms are minimap2's. MIT
license (see `LICENSE`), as for both projects.
