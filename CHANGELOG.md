# Changelog

Each release's binaries and full notes are on the
[releases page](https://github.com/DylanCosto/bonefish/releases). "Output" means the alignment
records. Each entry quotes the MD5 of the genome-wide HG001 sample's SAM output, `@PG` lines removed.

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
