# Changelog

Each release's binaries and full notes are on the
[releases page](https://github.com/DylanCosto/bonefish/releases). "Output" means the alignment
records. Each entry quotes the MD5 of the genome-wide HG001 sample's SAM output, `@PG` lines removed.

## Unreleased

**Fixed: base-modification tags on partial records.** With `-y`, a read's `MM`/`ML`/`MN` tags
were copied onto every record, including hard-clipped supplementaries and secondaries written
with SEQ `*`, where their positions (counted over the whole read) point past or into the wrong
bases. They are now kept only on records whose SEQ is the whole read (primaries, and
supplementaries with `-Y`). Output without `-y` is unchanged.

## v0.5.1 (2026-11-24)

**New: `--progress FILE`.** The progress the `--tui` panel shows, written to a file for a
pipeline's display (lrpipe's): one line, the fraction done and what is being done,
tab-separated (`0.4130	aligning · 2.4M reads · part 2/~4 · ~5m left`), replaced by a rename
as the run goes and ending at `1.0000	done · N reads`, with or without `--tui` and a terminal.
Without `--progress`, bonefish writes to the file `$LRPIPE_PROGRESS` names, as `lrpipe run` sets
it for its align stage. Output and logging are unchanged.

**Release check:** on 144,075 HG003 R10.4.1 reads (chr20) against GRCh38's chr20, SAM output
and sorted-BAM records are byte-identical to v0.5.0 (219,080 BAM records), with and without
`--progress`, at the same speed (69.9 s against 69.9 s for the sorted BAM, 16 threads).

## v0.5.0 (2026-11-18)

**New: live progress display (`--tui`).** One panel on stderr follows the run through its
stages (reference, alignment, index-part merge, sorted BAM, BAM index), each with its own clock:

```
bonefish  Running  00:29
┌─ Alignment progress ──────────────────────────────────────────────────────────────────┐
│  Prepare reference   done                                                       00:20 │
│› Align reads         [██████░░░░]  63% · 5252 reads · part 3/~4 · ~12s left     00:08 │
│  Sort and write BAM  queued                                                           │
│  Write BAM index     queued                                                           │
└───────────────────────────────────────────────────────────────────────────────────────┘
```

- **Alignment:** overall progress from how far into the read files bonefish is, reads aligned,
  reads/s and the time left. With a reference split into index parts (`-I`), each part is one
  pass over the reads: `part 3/4`, or `part 3/~4` while the reference is still being read. The
  current row stays on alignment while later parts are indexed.
- **Reference:** whether it is reading the FASTA, building the index or loading a prebuilt one.
- **Sorted BAM:** a progress bar while the BAM is written, then the record count.
- **Stays out of the way:** off unless asked for, and only when stderr is a terminal that
  alignment records are not printed to; otherwise, and with `TERM=dumb`, the usual status
  messages. Warnings print above the panel. Under 56 columns it is one line. Ctrl-C marks the
  run `Interrupted`.

**zflate replaces libdeflate.** BGZF input and output, the sorter's temporary runs and the BAI/CSI
indexes go through [zflate](https://github.com/DylanCosto/zflate), a pure-Rust port of libdeflate,
instead of libdeflate's C library. A build without jemalloc or mimalloc no longer needs a C
compiler. zflate needs Rust 1.91 or newer.

- **Output:** unchanged (`9bc7ff3e…`). Sorted BAMs, from FASTQ and from dorado BAM input, hold the
  same records, and their BAI and CSI indexes answer region queries as samtools' own do. At the
  default level the genome-wide sample's BAM is 465 bytes smaller.
- **Speed:** CPU time is unchanged.
- **`--sort-level`** outside 0-12 is rejected with an error, instead of a panic in a compression
  thread.

**Release check:** on 188,140 HG002 R10.4.1 reads (chr20) against GRCh38, SAM and sorted-BAM
output are byte-identical to v0.4.0, at the same speed (81.7 s against 81.6 s, 16 threads).
Clair3 and Sniffles2 calls on chr20 score the same against GIAB: SNP F1 0.9987, indel F1 0.9016,
SV F1 0.8139.

**Fixes**

- **Multi-part SAM/BAM output** (a reference larger than `-I`, 8G by default): the header listed
  only the first part's contigs, so alignments to the others pointed at missing contigs and BAM
  output failed. It now goes through the split-index merge, as `--split-prefix` does, mapped on all
  threads, with a complete header. Its temporary files go next to `--sort-tmp` or `-o`.
- **Multi-part SAM/BAM output reads the queries twice**, so it needs them as files: reads from
  stdin or a pipe are refused with an error (they used to be accepted with a warning, and with
  `--split-prefix` and a prebuilt index gave no records). Use a larger `-I`, or PAF output.
- **A reference that fits in one `-I` batch is one index part.** When its last sequence pushed the
  batch past `-I` (e.g. chr20 with `-I20M`), SAM/BAM output went through the split-index merge
  for that one part: slower (7.2 s against 6.0 s on 10,000 ONT reads), with MAPQ changed on about
  3% of records, and reads from stdin refused. The output now matches a run without `-I`.
- **Prebuilt indexes with a non-default k** (e.g. built with `-k15 -w10`) supply k for chaining as
  well as seeding, so they give the same output as mapping the FASTA with the same options.
- **Output errors end the run.** When the output could not be written (a closed pipe, a full
  disk), bonefish hung at 0% CPU if reads were left; it now exits with the error. SAM and PAF
  output is flushed before success is reported.
- **Ctrl-C removes temporary files.** On SIGINT, SIGTERM or SIGHUP bonefish removes the
  split-index directory and the sorter's temporary runs (`<prefix>.NNNN.run`) before it exits
  (exit status 128 + the signal); a sorted-BAM run interrupted mid-way used to leave hundreds of
  run files. On a FUSE filesystem such as ntfs-3g an empty split directory can remain.

## v0.4.0 (2026-10-28)

**New: `-x hifi`** (alias `map-hifi`), for PacBio HiFi genomic DNA. It is the `dna` preset, plus
one extra step: stretches of fragmented long gaps are realigned under minimap2 `map-hifi` scoring,
so one SV is no longer split into several false ones.

- **Speed:** on whole-genome HG002 it needs 3.9x less CPU time than minimap2 `map-hifi`.
- **Accuracy:** SNP and indel F1 are higher than minimap2's (0.9984 and 0.9922), and SV F1 is
  0.003 lower. Report: `docs/validation/2026-10-28-hifi/`.

**Preset names:** `ont` is now an alias for `dna`, and preset names ignore case. `dna` output
is unchanged on the golden-output tests.

## v0.3.0 (2026-10-21)

**New: `-x tx`.** The preset maps direct RNA or cDNA to a transcriptome for quantifiers such as
oarfish and NanoCount (minimap2 `-ax map-ont -N 100`).

- **Speed:** on HG002 it is 2.3x (cDNA) and 2.8x (direct RNA) faster than minimap2.
- **Accuracy:** it maps the same reads as minimap2. oarfish gene counts agree at r = 0.9995 or
  higher. Report: `docs/validation/2026-10-21-transcriptome/`.

**Faster with many reference sequences:** the ALT-contig count is computed once per index, not
by scanning every target for every read. This matters for transcriptomes and metagenome
assemblies.

**Faster chaining:** isolated anchors are dropped before chaining. DNA uses 1.2% less CPU, with
output unchanged (`9bc7ff3e…`).

## v0.2.0 (2026-10-21)

**New: direct RNA and cDNA.** Spliced alignment to a genome is available through two presets:

- `-x rna`: direct RNA, as minimap2 `-ax splice -uf -k14`.
- `-x cdna`: cDNA, as minimap2 `-ax splice`. `splice:hq` is also available.

Spliced runs accept `--junc-bed` and minimap2's splice options, and add the `ts:A` transcript-strand
tag.

**Speed:** against minimap2 on HG002, direct RNA runs 2.6x faster and cDNA 5.4x.

- Introns are aligned against only the reference near the chain's anchors, with a 10-mer search
  for exons that have no anchor.
- cDNA chains are aligned once instead of once per strand.
- Chaining uses the SIMD chainers, and read-end extensions use reduced windows.

**Accuracy:** 99.5-99.9% of primary alignments are identical to minimap2's, and GENCODE junction
accuracy is unchanged. Report: `docs/validation/2026-10-21-rna-cdna/`.

**DNA:** output unchanged (`9bc7ff3e…`).

**Removed:** the unused GPU experiment and the WebAssembly browser demo.

## v0.1.7 (2026-10-21)

**AVX-512 kernels** for the wide and narrow DP, chaining, small sorts and minimizer sketching,
chosen at run time.

- On a Xeon Gold 5418Y, the genome-wide sample takes 27% less CPU time than with v0.1.6.
- On AVX2 machines it takes 3% less.
- Output unchanged (`9bc7ff3e…`).

## v0.1.6 (2026-10-16)

**Fix:** the BAI/CSI index written with `-o x.bam` had about 2.5x too many chunks, so region readers
(Clair3, longphase, IGV) were slower on bonefish BAMs.

- The index is now identical in content to `samtools index`.
- Clair3 on HG002 chr20 runs in 181 s instead of 216 s.
- Output unchanged.

## v0.1.5 (2026-10-16)

**Sub-anchored gap filling is back**, outside tandem repeats, recovering most of v0.1.4's cost.

- Whole-genome re-validation on HG002 and HG005: SNP and indel F1 equal to minimap2's.
- The tandem-repeat benchmark matches minimap2.
- Output `9bc7ff3e9f498cc1530b2b1d78c47e10`.

## v0.1.4 (2026-10-15)

**Fix:** fewer false SV calls in tandem repeats.

- Sub-anchored gap filling is off: it could align reads through the wrong repeat copy.
- Costs about 6% CPU.
- Output `8fbcd8ce9df67b78a978ef6497873caa`.

## v0.1.3 (2026-10-15)

**Fixes:**

- Reads from a pipe now work.
- A CSI index is written for references longer than 512 Mbp.
- 7 crashes on corrupt BAM/BGZF input are fixed.

**Changes:** the default batch is `-K 100M`, using 2.7 GB less memory. Output unchanged.

## v0.1.2 (2026-10-14)

**Input and CLI hardening:**

- Malformed input is fatal, with a clear error.
- `-o x.bam` implies sorting.
- `-y` copies dorado BAM tags and `@RG` lines.

## v0.1.1 (2026-10-14)

**Platforms and validation:**

- Faster arm64 chaining.
- The DP kernel is ported to NEON and SSE4.1, so arm64 output equals x86.
- First whole-genome validation (HG002).
- Output `deed7949d35f15bd8eed935be47a1d81`.
