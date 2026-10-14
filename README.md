# bonefish

A fast aligner for Oxford Nanopore DNA reads, written in Rust. bonefish takes ONT reads (FASTQ,
FASTA, or the BAM that dorado writes) and produces SAM, PAF, or a coordinate-sorted, indexed BAM
in one step. It is tuned and validated for variant calling: on whole human genomes (GIAB HG002 and
HG005, R10.4.1), its SNP and indel accuracy match minimap2's, its structural-variant F1 is
slightly higher and its CpG methylation calls (modkit) agree with minimap2's at r = 0.999, while it goes from reads to a sorted, indexed BAM about 3x faster than
`minimap2 | samtools sort` (aligner alone: 2.6-8x).

## Install

Prebuilt binaries are attached to each GitHub release: Linux x86_64 and arm64 (glibc 2.28 or
newer: RHEL/Rocky 8+, Ubuntu 20.04+, Debian 10+) and macOS on Apple silicon. Unpack and run
`./bonefish`. To build from source with a recent stable Rust toolchain:

```bash
cargo build --release --features jemalloc    # jemalloc is ~10% faster than the system allocator
```

The SIMD kernels are chosen at run time (AVX-512/AVX2 on x86_64, NEON on arm64), so one binary
runs on any CPU of its architecture. The speed and accuracy figures below were measured on x86_64
with AVX2. The release workflow builds each platform's binary only after that platform passes
the golden-output regression tests.

## Quick start

```bash
bonefish -d ref.raw --idx-raw ref.fa                  # optional: build a fast-loading index
bonefish -t 16 -y -o sample.bam ref.raw reads.bam     # dorado BAM in, sorted BAM + .bai out
bonefish -t 16 -a ref.fa reads.fq.gz > sample.sam     # FASTQ in, SAM out
```

- There is one preset, `dna` (ONT genomic DNA, R10.4.1 or R9.4.1), applied by default. The
  minimap2 names `lr:hq`, `lr:hq-fast`, `map-ont` and `lr` are accepted for it; other minimap2
  presets (PacBio, short reads, RNA/splice, assembly) are not supported.
- Reads can be FASTQ, FASTA or BAM, plain or gzipped, from files or stdin; several files are
  mapped in turn. BAM input may be unaligned (dorado) or aligned: reverse-strand records are
  restored to the read's orientation and secondary/supplementary records are skipped.
- `-y` copies each read's tags to the output: from a dorado BAM, everything except old alignment
  tags (base modifications `MM`/`ML`/`MN`, `RG`, `qs`, `ch`, `st`, ...), plus the input's `@RG`
  header lines (which carry the basecalling model); from FASTQ, the SAM-style tags of the comment.
  `-R` sets the read group instead.
- `-o FILE.bam` writes a coordinate-sorted BAM and `FILE.bam.bai` directly (same records and order
  as `bonefish -a | samtools sort --write-index`); `--sort-mem` (default 4G) sets the sorter's
  in-memory batches, with spills to `--sort-tmp`. A human genome needs about 15-20 GB with 16
  threads, and about 23 GB when `-y` carries methylation tags into a sorted BAM; `--sort-mem 1G`
  brings that to about 19 GB on smaller machines.
- Malformed input (truncated records, quality/sequence length mismatches, read names over the
  SAM limit, non-sequence files) stops the run with an error naming the file and line, and a
  non-zero exit status.
- Most minimap2 command-line options (`-t`, `-a`, `-c`, `--cs`, `--MD`, `--eqx`, `-N`, `-p`, ...)
  work as in minimap2.

## Speed and accuracy (16 threads)

Whole genomes, reads to sorted and indexed BAM (minimap2 2.31 `-ax lr:hq | samtools sort`):

| sample (R10.4.1 sup) | minimap2 + samtools | bonefish | SNP F1 | indel F1 | SV F1 (v5.0q) |
|---|---:|---:|---|---|---|
| HG005, 27.9x | 114.7 min | **34.8 min** | 0.9975 / 0.9974 | 0.9395 / 0.9396 | – |
| HG002, 27.2x | 114.8 min | **37.5 min** | 0.9971 / 0.9971 | 0.8737 / 0.8737 | 0.8055 / **0.8094** |

(F1 pairs are minimap2 / bonefish; Clair3 vs GIAB v4.2.1, Sniffles2 vs GIAB v5.0q.) Aligner alone:

| dataset | minimap2 `lr:hq` | bonefish | speedup |
|---|---:|---:|---:|
| HG002 R10.4.1 HAC, chr20 reads (3.5 Gbp) | 417 s | 72 s | 5.8x |
| HG002 R9.4.1, chr20 reads (4.6 Gbp) | 412 s | 103 s | 4.0x |
| HG002 R10.4.1 ultra-long (3.4 Gbp) | 542 s | 65 s | 8.3x |
| HG001 R10.4.1 sup, random genome-wide sample (2.7 Gbp) | 202 s | 44 s | 4.6x |
| *Salmonella enterica* R10.4.1 (1.2 Gbp) | 25 s | 10 s | 2.6x |

Chromosome-scale variant calling (Clair3 + rtg vcfeval vs GIAB v4.2.1; Sniffles2 + Truvari vs GIAB v5.0q):

| dataset | aligner | SNP F1 | indel F1 | SV F1 |
|---|---|---:|---:|---:|
| HG002 HAC chr20 | minimap2 | 0.9985 | 0.8294 | 0.7887 |
| | bonefish | 0.9984 | 0.8310 | 0.7917 |
| HG002 R9.4.1 chr20 | minimap2 | 0.9980 | 0.8173 | 0.7704 |
| | bonefish | 0.9980 | 0.8172 | 0.7734 |

Against other ONT aligners on the genome-wide sample (same machine, 16 threads): dorado 2.1.2
aligner 4.3x slower to SAM and 10.5x slower to a sorted, indexed BAM; mm2plus 5.1x, mm2-fast 14x
and Winnowmap2 38x slower; NGMLR about 260x slower on a 10% subset. Datasets, commands, memory,
the full comparison and the known differences from minimap2 are in
[docs/benchmarks.md](docs/benchmarks.md); the raw measurements behind them are in
[docs/validation/](docs/validation/) (whole genomes:
[2026-10-12-giab-wgs](docs/validation/2026-10-12-giab-wgs/README.md); other aligners:
[2026-10-14-aligner-comparison](docs/validation/2026-10-14-aligner-comparison/README.md);
methylation: [2026-10-14-modkit-methylation](docs/validation/2026-10-14-modkit-methylation/README.md)).

**Apple Silicon.** arm64 builds produce the same output as x86: the genome-wide HG001 sample gives
the identical MD5 (`deed7949d35f15bd8eed935be47a1d81`) on an Apple M4 Pro and on x86
([report](docs/validation/2026-10-13-hg001-macos-parity/README.md)). The NEON chaining and DP
kernels and their measurements on HG002 ultra-long data are documented in
[the macOS report](docs/validation/2026-10-13-macos-arm64/README.md) and
[the clean-DP report](docs/validation/2026-10-13-neon-clean-dp/README.md), with
[reproduction scripts](scripts/validation/README.md).

## How it differs from minimap2

Same seeding and scoring as minimap2's `lr:hq`, with cheaper alignment and chaining: narrower,
adaptive DP bands with retries where an indel needs a wider band, a probe that stops extending
into low-identity read ends, ungapped-first and sub-anchored gap filling, a 5 kb reference-gap limit
in the main chaining pass (with a long-join rescue that skips repeat-dense reads), segment-tree
RMQ chaining and diagonal-bucketed exact chaining, plus multithreaded I/O and a built-in sorter.
Output is therefore not byte-identical to minimap2. `bonefish --help` lists the individual
options.

## Limitations

- DNA only: no spliced (RNA/cDNA) alignment, and no PacBio, short-read or assembly presets.
- Output is SAM, PAF or sorted BAM (no CRAM).
- Indexes written with `-d` are bonefish's own format; minimap2 `.mmi` files can be read but not
  written.
- A human genome needs about 20 GB of memory with 16 threads (about 23 GB for a sorted BAM with
  methylation tags; see `--sort-mem`).
- The browser demo (`web/`, WebAssembly) still uses the older DP kernels, so its alignments can
  differ slightly from the native builds.

## Library

`bonefish-core` (imported as `bonefish`) exposes the aligner as a Rust library
(`bonefish::Aligner`, `Preset::Dna`); see `bonefish-core/examples/`.

## Credits and license

bonefish is a fork of [rammap](https://github.com/jwanglab/rammap), a pure-Rust reimplementation
of [minimap2](https://github.com/lh3/minimap2) by Heng Li. The index format, SIMD DP kernels and
most of the mapping pipeline come from rammap; the alignment algorithms are minimap2's. MIT
license (see `LICENSE`), as for both projects.
