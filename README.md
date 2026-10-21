# bonefish

A fast aligner for Oxford Nanopore reads, written in Rust. bonefish maps genomic DNA, direct
RNA and cDNA reads to a reference genome. It writes SAM, PAF or a sorted, indexed BAM in one step.

- **Accurate.** On whole human genomes (GIAB HG002 and HG005), SNP, indel and structural-variant
  calls match minimap2's. Methylation calls agree at r = 0.999.
- **Fast.**
  - Aligner alone: 2.6-8x faster than minimap2 on DNA, 2.6x on direct RNA and 5.4x on cDNA.
  - Reads to a sorted BAM: about 3x faster than `minimap2 | samtools sort`.
- **Drop-in.** It reads dorado's BAM output, keeping methylation and other tags, and accepts
  minimap2's options and preset names.

## Install

**Prebuilt binaries:** download from the [latest release](https://github.com/DylanCosto/bonefish/releases/latest).
They are available for Linux x86_64 and arm64 (glibc 2.28+: RHEL/Rocky 8+, Ubuntu 20.04+, Debian
10+) and macOS on Apple silicon. Unpack and run `./bonefish`.

**From source**, with a recent stable Rust toolchain:

```bash
cargo build --release --features jemalloc
```

**Container:**

```bash
docker build -t bonefish https://github.com/DylanCosto/bonefish.git
```

## Quick start

```bash
# DNA: dorado BAM in, sorted BAM + index out (-y keeps methylation tags)
bonefish -t 16 -y -o sample.bam ref.fa reads.bam

# Direct RNA and cDNA, spliced to the genome
bonefish -t 16 -x rna  -y -o rna.bam  ref.fa reads.bam
bonefish -t 16 -x cdna -y -o cdna.bam ref.fa reads.bam

# FASTQ in, SAM out
bonefish -t 16 -a ref.fa reads.fq.gz > sample.sam
```

For repeated runs, build an index once: `bonefish -d ref.raw --idx-raw ref.fa`. It loads in under
a second; add `-x rna` or `-x cdna` for a spliced index. New to bonefish? The
[tutorial](docs/tutorial.md) goes from ONT's public demo data to a sorted BAM and methylation calls
in a few minutes.

| preset | for |
|---|---|
| `dna` (default) | genomic DNA: R10.4.1 or R9.4.1, sup or hac, ultra-long (minimap2 `lr:hq`) |
| `rna` | direct RNA (minimap2 `splice -uf -k14`) |
| `cdna` | cDNA (minimap2 `splice`) |

Input formats, output options, memory use and CPU support are covered in
[docs/usage.md](docs/usage.md).

## Speed and accuracy

All runs below used 16 threads and compare against minimap2 2.31 on the same reads.

**Whole human genomes**, reads to a sorted, indexed BAM (v0.1.5; HG005 timed with v0.1.1). Each F1 cell gives minimap2 / bonefish:

| sample (R10.4.1 sup) | minimap2 + samtools | bonefish | SNP F1 | indel F1 | SV F1 |
|---|---:|---:|---|---|---|
| HG002, 27x | 115 min | **38 min** | 0.9971 / 0.9971 | 0.8737 / 0.8737 | 0.8055 / 0.8053 |
| HG005, 28x | 115 min | **35 min** | 0.9975 / 0.9975 | 0.9395 / 0.9395 | – |

Small variants were called with Clair3 against GIAB v4.2.1, and structural variants with
Sniffles2 against GIAB v5.0q.

**Aligner alone**:

| dataset | minimap2 | bonefish | speedup |
|---|---:|---:|---:|
| DNA, HG002 ultra-long (3.4 Gbp) | 542 s | 65 s | 8.3x |
| DNA, HG001 genome-wide sample (2.7 Gbp) | 202 s | 44 s | 4.6x |
| DNA, HG002 R9.4.1 chr20 (4.6 Gbp) | 412 s | 103 s | 4.0x |
| cDNA, HG002 PCS114 (809k reads) | 1,223 s | 227 s | 5.4x |
| direct RNA, HG002 RNA004 (902k reads) | 456 s | 177 s | 2.6x |

On RNA, 99.5-99.9% of primary alignments are identical to minimap2's, and junction accuracy
against GENCODE is the same.

bonefish is also faster than dorado's aligner (4.3x), mm2plus (5.1x), mm2-fast (14x) and
Winnowmap2 (38x).

- **Datasets and comparisons:** [docs/benchmarks.md](docs/benchmarks.md) has the datasets,
  commands, every comparison, and the known differences from minimap2.
- **Raw evidence:** [docs/validation](docs/validation/README.md) indexes the reports behind each
  number.
- **Other read types:** the reports cover metagenomes, duplex reads, other references
  (T2T-CHM13, ALT contigs), simulated reads with known origin, medically relevant genes, tandem
  repeats and output checks.

## How it works

bonefish uses minimap2's seed-chain-align method with the same scoring, so its output stays
close to minimap2's. Most of the speed comes from doing less dynamic-programming (DP) alignment:

- **DP bands:** narrower, adaptive bands, widened only where an indel needs it.
- **Gap filling:** gaps are first tried without indels, and long gaps are split at unique
  matches.
- **Introns:** they are aligned against just the reference next to the chain's anchors. A
  10-mer search looks for exons that have no anchor.
- **cDNA strands:** each chain is aligned once for both transcript strands.

The rest comes from engineering:

- **SIMD kernels:** hand-written AVX-512, AVX2 and NEON kernels for chaining, DP, seeding and
  sorting, chosen at run time.
- **I/O and sorting:** multithreaded I/O, an index that loads in under a second, and a built-in
  BAM sorter.

## Limitations

- ONT reads only: there are no PacBio, short-read or assembly presets.
- The RNA presets are new (v0.2.0). They are validated against minimap2, GENCODE and simulated
  reads, but not yet with downstream isoform or expression tools.
- Output is SAM, PAF or sorted BAM; there is no CRAM output.
- On low-accuracy reads (about 93% identity), minimap2's `map-ont` maps about 1% more reads.

## Library

`bonefish-core` (imported as `bonefish`) exposes the aligner as a Rust library: `bonefish::Aligner`
with `Preset::Dna`, `Preset::Rna` and `Preset::Cdna`. See `bonefish-core/examples/`.

## Credits and license

bonefish is a fork of [rammap](https://github.com/jwanglab/rammap), a pure-Rust reimplementation
of [minimap2](https://github.com/lh3/minimap2) by Heng Li. The index format, SIMD DP kernels and
most of the mapping pipeline come from rammap, and the alignment algorithms are minimap2's.
bonefish is MIT licensed (see `LICENSE`), as both of those projects are.
