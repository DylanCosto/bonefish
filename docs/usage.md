# Using bonefish

## Presets

| preset | for | minimap2 equivalent |
|---|---|---|
| `dna` (default; alias `ont`) | ONT genomic DNA: R10.4.1 or R9.4.1, sup or hac, ultra-long | `-x lr:hq`, with cheaper chaining and alignment |
| `hifi` | PacBio HiFi genomic DNA; uses `dna` indexes | `-x map-hifi` |
| `rna` | direct RNA (RNA004), spliced to a genome | `-x splice -uf -k14` |
| `cdna` | cDNA, spliced to a genome, transcript strand from splice signals | `-x splice` |
| `splice:hq` | accurate cDNA | `-x splice:hq` |
| `tx` | direct RNA or cDNA to a transcriptome, for quantification (oarfish, NanoCount) | `-x map-ont -N 100` |

The minimap2 names `lr:hq`, `lr:hq-fast`, `map-ont` and `lr` are accepted for `dna`, `map-hifi` for
`hifi` and `splice` for `cdna`. Preset names ignore case. PacBio Iso-Seq/Kinnex RNA has not been
validated. Short-read and assembly presets are not supported.

## Indexes

bonefish can map straight from a FASTA file, building the index in memory first. For repeated runs,
build an index once:

```bash
bonefish -d ref.raw --idx-raw ref.fa              # dna; loads in under a second (8.5 GB for GRCh38)
bonefish -x cdna -d cdna.raw --idx-raw ref.fa     # spliced presets need their own index (k/w differ)
bonefish -x rna  -d rna.raw  --idx-raw ref.fa
bonefish -x tx   -d tx.raw   --idx-raw transcripts.fa
```

Without `--idx-raw`, `-d` writes a smaller index that takes a few seconds to load. Indexes are
bonefish's own format. minimap2 `.mmi` files can be read but not written.

## Input

- **Formats:** FASTQ, FASTA or BAM, plain or gzipped, from files, stdin (`-`) or a pipe. Several
  files are mapped in turn.
- **BAM input:** may be unaligned (as dorado writes it) or aligned. Reverse-strand records are
  restored to the read's orientation, and secondary and supplementary records are skipped.
- **RNA:** `U` in direct-RNA reads is read as `T`.
- **Malformed input stops the run** with an error naming the file and a non-zero exit status. This
  covers truncated records, quality/sequence length mismatches, read names over the SAM limit,
  corrupt gzip, BGZF or BAM data, and non-sequence files.

## Output

- **SAM** with `-a`, and **PAF** by default. With `-c`, PAF carries CIGAR strings; `--cs`, `--MD`
  and `--eqx` work as in minimap2.
- **`-o FILE.bam`** writes a coordinate-sorted BAM and its `FILE.bam.bai` directly:
  - the same records and order as `bonefish -a | samtools sort --write-index`;
  - a `.csi` index instead when a reference is longer than 512 Mbp, as many plant chromosomes are;
  - `--sort-mem` (default 4G) sets the sorter's in-memory batches, which spill to `--sort-tmp`.
- **`-y` copies each read's tags:**
  - from a dorado BAM, everything except old alignment tags: base modifications `MM`/`ML`/`MN`,
    `RG`, `qs`, `ch`, `st` and so on, plus the input's `@RG` header lines, which carry the
    basecalling model;
  - from FASTQ, the SAM-style tags of the comment.

  `-R` sets the read group instead.
- **Spliced presets:**
  - add the `ts:A` transcript-strand tag;
  - `--write-junc` writes the junctions found;
  - `--junc-bed FILE` (BED12 transcripts, or one intron per BED6 line) gives annotated junctions a
    bonus.

  `-u`, `-G` (maximum intron length), `-C`, `--junc-bonus`, `--junc-pen` and `--splice-flank`
  work as in minimap2.

Most other minimap2 options (`-t`, `-N`, `-p`, `-K`, ...) work as in minimap2. `bonefish --help`
lists them all.

## Memory

With 16 threads, mapping a human genome (DNA) peaks at:

| output | peak memory |
|---|---:|
| SAM | about 13 GB |
| sorted BAM | 15-20 GB |
| sorted BAM, `-y` carrying methylation tags | about 23 GB |
| sorted BAM, `-y` with methylation tags, `--sort-mem 1G` | about 19 GB |

On a 16 GB machine, `-t 4 -K 10M` maps with a peak under 9 GB and gives the same output
([measurements](validation/2026-10-14-scaling-memory/README.md)).

The spliced presets' indexes are larger (11-13 GB for GRCh38): mapping RNA or cDNA with 16
threads peaks at 13-16 GB.

## CPUs

- **Runtime kernel choice:** the SIMD kernels are chosen at run time (AVX-512, AVX2 or SSE4.1 on
  x86_64; NEON on arm64), so one binary runs on any CPU of its architecture.
- **Same output everywhere:** the output is identical on every code path and across thread
  counts. arm64 gives the same output as x86
  ([Apple M4 Pro report](validation/2026-10-13-hg001-macos-parity/README.md)).
- **AVX-512 is faster:** on AVX-512 CPUs the DNA preset needs about 27% less CPU time than with
  AVX2 ([details](validation/2026-10-17-avx512/README.md)).
- **Forcing a code path** (testing only): `BONEFISH_FORCE_AVX2=1`, `BONEFISH_FORCE_SSE=1` and
  `BONEFISH_FORCE_SCALAR=1`.
