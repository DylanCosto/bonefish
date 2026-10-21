# Transcriptome mapping for quantification (21 Oct 2026)

Long-read quantifiers such as oarfish and NanoCount take reads aligned to a transcriptome, with
many secondary alignments per read. Their usual input is `minimap2 -ax map-ont -N 100`. This
report adds the `-x tx` preset for that workflow, and compares it with minimap2 and with
bonefish's `dna` preset.

## Data and method

- **Reads:** the first 200,000 HG002 direct-RNA (RNA004 sup) and cDNA (PCS114 hac) reads of the
  [RNA report](../2026-10-21-rna-cdna/README.md).
- **Simulated reads:** Badread reads from 30,000 GENCODE transcripts, cDNA-like and
  direct-RNA-like. Each read's source transcript and gene are known.
- **Reference:** GENCODE v49 transcripts, 533,740 sequences
  (`gencode.v49.transcripts.fa.gz`, sha256 `08ee7ee8…`).
- **Commands:**
  - minimap2 2.31: `-ax map-ont -N 100 --eqx`.
  - bonefish: `-a --eqx`, with either `-x tx` or `-N 100` on the `dna` preset.
  - Quantification: oarfish 0.10.3 with `--filter-group nanocount-filters`, compared with
    `tx_quant_compare.py`.
- **Machine:** Core i5-12600K, 16 threads.

**oarfish and bonefish BAMs:** oarfish only reads BAMs whose `@PG` header names an aligner it
has validated (minimap2, pbmm2, bramble or rammap). bonefish's BAMs were relabelled as rammap,
bonefish's upstream, for this test. oarfish otherwise only needs the standard `AS` tag.

## The `tx` preset

- **Seeding:** minimap2's `map-ont` seeding. k=15 and w=10, with no occurrence range and a 5 kb
  query gap limit, suit short, noisy RNA reads.
- **Secondary alignments:** up to 100 per read (`-N 100`).
- **Everything else:** the `dna` preset's chaining and DP. There are no introns in transcriptome
  mapping, so they apply unchanged.

## Results

| 200k reads | minimap2 `map-ont -N 100` | bonefish `dna -N 100` | **bonefish `-x tx`** |
|---|---:|---:|---:|
| cDNA: wall / CPU-s | 185 s / 2,285 | 146 s / 1,535 | **81 s / 1,149** |
| direct RNA: wall / CPU-s | 106 s / 1,026 | 102 s / 1,032 | **37 s / 485** |
| unmapped reads (cDNA / direct RNA) | 1,987 / 37,838 | 2,493 / 39,115 | **1,987 / 37,838** |
| oarfish genes vs minimap2: Pearson (cDNA / direct RNA) | | 0.99937 / 0.99968 | 0.99950 / 0.99997 |
| reads assigned to another gene (cDNA / direct RNA) | | 0.65% / 0.55% | 0.66% / 0.25% |
| reads assigned to another transcript (cDNA / direct RNA) | | 1.63% / 1.27% | 1.74% / 0.93% |
| simulated cDNA: primary on true transcript / gene | 93.73% / 99.53% | 93.89% / 99.48% | 93.88% / 99.53% |
| simulated direct RNA: mapped | 99.69% | 98.82% | 99.69% |
| simulated direct RNA: primary on true transcript / gene | 91.98% / 99.33% | 91.33% / 98.49% | 92.17% / 99.34% |
| simulated: true transcript among alignments (cDNA / direct RNA) | 99.89% / 99.68% | 99.85% / 98.80% | 99.89% / 99.68% |

- **Sensitivity:** `-x tx` leaves exactly the same number of reads unmapped as minimap2. The `dna`
  preset's sparser k=19 seeds missed some short reads.
- **Accuracy:** on simulated reads, `-x tx` puts at least as many primaries on the true
  transcript as minimap2. The true transcript is among the alignments as often.
- **Quantification:** at gene level, quantification from `-x tx` agrees with minimap2's at
  r = 0.9995 (cDNA) and 0.99997 (direct RNA). Assignment between isoforms, which is close to a
  tie for many reads, differs for 0.9-1.7% of reads.
- **Speed:** `-x tx` is 2.3x (cDNA) and 2.8x (direct RNA) faster than minimap2 in wall time.

## What made it fast

The first `tx` build was no faster than minimap2 (`runs/tx1`: cDNA 184 s, direct RNA 102 s).
Profiling found three things:

1. **A per-read scan of every reference sequence.** The code counted ALT contigs by looking at
   all 533,740 transcript names twice per read. That was 32 of 35 post-chain thread-seconds on
   20k reads. On the 195-sequence human genome it cost nothing, so it had never shown up. The
   count is now made once per index (`runs/tx3`: cDNA 85 s, direct RNA 39 s). Anything with many
   reference sequences, such as metagenome assemblies, gains as well.
2. **Isolated anchors.** A read's k-mers hit hundreds of short transcripts, and most hits have no
   other hit on the same transcript. They are now dropped before chaining, because they can only
   form one-anchor chains, which are always discarded. Output is identical. Chaining is 8%
   faster here, and DNA uses 1.2% less CPU (`runs/tx4`).
3. **Anchor copies.** Chain anchors are now copied only for regions that survive the secondary
   filter. This made little difference on its own.

For all three changes, output stayed identical: the DNA genome-wide sample has the same MD5
(`9bc7ff3e…`), every golden-output test passes, and the transcriptome SAM is unchanged.

## Files

- `evidence/runs/<run>/`: GNU time, stderr, oarfish logs and metadata.
  - `mm2`: minimap2.
  - `dna_N100`: the `dna` preset with `-N 100`.
  - `tx1`, `tx3`, `tx4`: the `tx` preset at each step above.
- `evidence/ledger.tsv`: one line per run.
- `evidence/scripts/`: `tx_bench.sh` (runs and scoring), `tx_quant_compare.py` (oarfish count
  correlation) and `tx_sim_truth.py` (simulated-read truth).

The simulated reads are those of the RNA report (`rna_sim_reads.sh`). `SHA256SUMS` covers
`evidence/`. Local paths are replaced by `${BENCH}`, `${DATA}`, `${REPO}` and `${HOME}`.
