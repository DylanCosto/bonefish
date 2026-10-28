# PacBio HiFi genomic DNA (28 Oct 2026)

This report adds the `-x hifi` preset and compares it with minimap2 2.31 `map-hifi`.

## Data and method

- **Reads:** PacBio Revio 2022Q4 HiFi data.
  - Whole genome: HG002 rep1, 6.53M reads, about 30x.
  - chr20: reads with a primary alignment on chr20 in PacBio's own BAMs, for HG002, HG003 and
    HG004. These reads are aligned to the whole genome.
- **Reference:** GRCh38 no-alt analysis set.
- **Small variants:** Clair3 (`hifi_revio` model), scored with `rtg vcfeval` against GIAB v4.2.1.
- **Structural variants (HG002):** Sniffles2, scored with Truvari (`--pctseq 0`) against GIAB v5.0q
  stvar, with strict and `--pick multi` matching.
- **Commands:**
  - minimap2: `-ax map-hifi` (`lr:hq` also run on chr20).
  - bonefish: `-a`, with either the `dna` preset or `-x hifi`.
  - All runs: `-y -t 16`, output piped to `samtools sort`.
- **Machine:** Core i5-12600K, 16 threads.

## The `hifi` preset

The `dna` preset scores alignments like minimap2's `lr:hq`: match 2, long-gap extension 1. Under
that scoring, a long deletion costs about as much as aligning through it on a few chance matches.
The deletion is then reported as several deletions with short aligned stretches between them.
On HiFi reads nearly every long gap is a real SV, so one deletion becomes several false SV calls.

`-x hifi` is the `dna` preset plus one step at output:

1. **Find the stretches:** two or more gaps of at least 50 bp, with no solid run between them. A
   solid run is a match of at least 50 bp at 90% identity or more.
2. **Realign them:** each stretch is realigned globally under `map-hifi` scoring
   (`-A1 -B4 -O6,26 -E2,1`).
3. **Merge what is left:** if the result is still split and the middle is below 90% identity,
   the middle becomes one deletion plus one insertion.

The step is skipped for stretches in tandem repeats and at alignment ends. Placement, MAPQ and
filtering are computed before it, so every read lands where `dna` puts it, with the same MAPQ.
`-x hifi` uses `dna` indexes.

## Results

**Whole genome, HG002:**

| | minimap2 `map-hifi` | bonefish `dna` | **bonefish `-x hifi`** |
|---|---:|---:|---:|
| aligner CPU time | 91,927 s | 23,110 s | **23,633 s** |
| reads to sorted BAM, wall time | 7,908 s | 3,438 s | **3,010 s** |
| peak memory | 18.7 GB | 14.5 GB | 14.3 GB |
| SNP F1 (FP) | 0.9980 (8,163) | 0.9984 (5,818) | **0.9984 (5,705)** |
| indel F1 (FP) | 0.9868 (5,147) | 0.9922 (3,694) | **0.9922 (3,642)** |
| SV F1, strict / multi | **0.8091 / 0.8748** | 0.8026 / 0.8670 | 0.8061 / 0.8704 |

**chr20** (SNP F1 / indel F1; SVs strict / multi for HG002):

| sample | minimap2 `map-hifi` | minimap2 `lr:hq` | bonefish `dna` | **bonefish `-x hifi`** |
|---|---|---|---|---|
| HG002 | 0.9971 / 0.9837 | 0.9974 / 0.9884 | 0.9975 / 0.9885 | **0.9975 / 0.9885** |
| HG002 SVs | 0.8110 / 0.8923 | 0.8082 / 0.8888 | 0.8076 / 0.8889 | 0.8096 / 0.8910 |
| HG003 | 0.9981 / 0.9803 | 0.9983 / 0.9836 | 0.9989 / 0.9830 | **0.9989 / 0.9830** |
| HG004 | 0.9981 / 0.9784 | 0.9984 / 0.9816 | 0.9984 / 0.9807 | **0.9986 / 0.9813** |

- **Small variants:** on every sample, `-x hifi` is better than `map-hifi` on both SNPs and
  indels, and equal to or better than `dna`.
- **SVs:** `-x hifi` recovers about half of the gap between `dna` and `map-hifi`. The remaining
  0.003 F1 (strict) mostly comes from reads that minimap2 splits into two alignments, where
  bonefish keeps one.
- **Speed:** `-x hifi` needs 3.9x less CPU than `map-hifi`, and gets reads to a sorted BAM 2.6x
  faster. The realignment costs about 2% CPU over `dna`.

**ONT reads:** the same realignment was tested on HG002 R10.4.1 sup reads (chr20). It added 2 true
SVs (strict F1 0.8139 to 0.8156, equal to minimap2 `lr:hq`), and small-variant calls did not
change. `dna` is unchanged in this release, and `ont` is now an alias for it.

## Files

- `evidence/<set>/<aligner>/`: `accuracy.tsv` (vcfeval), `sv_accuracy.tsv` (Truvari),
  `align.time` (GNU time: wall, user, system, max RSS in KB) and `compare_vs_mm2hifi.txt` (primary
  placement against minimap2 `map-hifi`). The sets are `wgs_HG002`, `chr20_HG00{2,3,4}` and
  `ont_chr20_HG002`.
- `evidence/scripts/`: data download (`hifi_fetch.sh`, `hifi_wgs_fetch.sh`) and runs
  (`hifi_eval.sh`, `hifi_wgs_eval.sh`, `ont_realign_eval.sh`).

`SHA256SUMS` covers `evidence/`. Local paths are replaced by `${BENCH}`, `${DATA}` and `${HOME}`.
