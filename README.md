# bonefish

A fast aligner for Oxford Nanopore DNA reads, written in Rust. bonefish takes ONT reads (FASTQ,
or the unaligned BAM that dorado writes) and produces SAM, PAF, or a coordinate-sorted, indexed BAM
in one step. It is tuned and validated for variant calling: against GIAB truth sets, its SNP and
indel accuracy match minimap2's and its structural-variant recall is at or above minimap2's, while
aligning 2.6-8x faster.

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

Apple Silicon chaining experiments and the nominal 1×/10× HG002 ONT WGS runs are documented in
[the macOS validation report](docs/validation/2026-10-13-macos-arm64/README.md), with source
provenance, trial logs, output checks and [reproduction scripts](scripts/validation/README.md).
The Mac output-preservation checks do not establish real-data x86/ARM equality or new accuracy
qualification; those limits are recorded with the measurements.

The [clean NEON/SSE4.1 DP report](docs/validation/2026-10-13-neon-clean-dp/README.md)
records native platform checks and a paired nominal 1× HG002 experiment: 16.85% lower
mean elapsed time on the tested Mac, with all six trial receipts and changed-read evidence.
Older ARM alignment records change; full WGS x86/ARM parity remains untested.

A subsequent [HG001 genome-wide Mac replay](docs/validation/2026-10-13-hg001-macos-parity/README.md)
reproduced the supplied x86 MD5 prefix with 161,762 reads. The report records the full
Mac MD5/SHA256 and structural validation; full-value x86 confirmation remains pending.

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
