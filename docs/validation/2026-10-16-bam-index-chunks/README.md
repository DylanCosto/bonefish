# Why Clair3 ran slower on bonefish BAMs (16 Oct 2026)

Clair3 took 9-15% longer on bonefish's BAMs than on minimap2's in the whole-genome runs
([12 Oct](../2026-10-12-giab-wgs/README.md), [15 Oct](../2026-10-15-giab-wgs-v2/README.md)). The
cause was the BAI index that `bonefish -o x.bam` writes, not the alignments.

## Where the time went

Clair3's job logs (`evidence/clair3_job_logs/`), HG002 whole genome, wall time per stage:

| Clair3 stage | minimap2 | bonefish v0.1.1 | bonefish v0.1.5 |
|---|---:|---:|---:|
| 1. pileup calling | 1,564 s | 1,803 s | 2,009 s |
| 3. phasing (longphase) | 99 s | 416 s | 434 s |
| 6. full-alignment calling | 2,946 s | 2,940 s | 2,912 s |

HG005 shows the same pattern (phasing 99 s against 423 s). Both slow stages read the BAM by region
through its index; full-alignment calling takes the same time.

## Cause

On chr20 (HG002 PAW70337, 55x), the same bonefish records were written three ways and phased with
longphase as Clair3 runs it (`evidence/timings.txt`):

| BAM | longphase |
|---|---:|
| bonefish v0.1.5 (`-o x.bam`, its own sorter and BAI) | 46.1 s (142 CPU-s) |
| the same BAM, re-indexed with `samtools index` | 10.5 s |
| the same records through `samtools sort --write-index` | 10.4 s |
| minimap2 2.31, `samtools sort` | 10.7 s |

The BAM itself was fine; the index split each bin into far more chunks than htslib's: 168,634
chunks against 66,877 for chr20 (`evidence/index_comparison.txt`). A region query reads chunk by
chunk, so it sought to and inflated the same 64 KB block once per small chunk. Two differences from
htslib's indexer were responsible:

1. **Chunk ends at block boundaries.** A record that ends exactly where its BGZF block's data ends
   was given the end offset "block N, end of data", where readers and htslib give "block N+1, offset
   0" (htslib's `hts_idx_amend_last`). The next record of the same bin, starting at block N+1, then
   never looked contiguous, so every block boundary started a new chunk.
2. **No `compress_binning`.** htslib finishes an index by moving bins whose chunks span less than
   64 KB of compressed file into their parent bin, and merging chunks that start in the block where
   the previous one ends.

A third, minor difference: linear-index windows that no record overlaps got the previous window's
offset; htslib gives them the next record's, which is tighter.

## Fix

The sorter's index writer (`bonefish-core/src/sorted_bam.rs`) now does all three as htslib does.
Its BAI is identical in content to `samtools index` of the same BAM (every bin, chunk and
linear-index entry; only the order bins are stored in differs, which readers ignore), and so is
its CSI to `samtools index -c` (checked on a 664 Mbp contig, past BAI's limit). The BAM is
unchanged.

| chr20, `run_clair3.sh --ctg_name chr20`, 16 threads | wall | CPU-s |
|---|---:|---:|
| bonefish v0.1.5 index | 216 s | 1,252 |
| **fixed index** | **181-182 s** | 1,165-1,180 |
| samtools index of the same BAM | 177 s | 1,142 |
| minimap2 2.31 BAM | 174-179 s | 1,148-1,210 |

The calls from the three bonefish indexes are byte-identical (145,198 records). With the fixed index,
Clair3 on bonefish's BAM runs within run-to-run variation of minimap2's, and longphase alone is as
fast (10.6 s against 10.7 s). The whole genome was not re-run; at 4.5x faster phasing and the pileup
stage's region reads no longer repeated, the 9-15% whole-genome gap should close the same way.

BAMs written by v0.1.5 and earlier are valid; re-indexing them with `samtools index` gives the
faster index.

## Tests

`sorted_bam_output` (regression test) now checks that no chunk in the written BAI ends at the end
of a block's data; it fails on the v0.1.5 writer. Unit tests cover the bin and chunk merging and the
linear-index fill against htslib's behaviour.

## Files

`evidence/timings.txt`, `evidence/index_comparison.txt`, `evidence/compare_index.py` (compares two
BAI or CSI files by content), `evidence/clair3_job_logs/` (Clair3's GNU parallel job logs for the
pileup, phasing and full-alignment stages of the whole-genome runs). `SHA256SUMS` covers
`evidence/`. Local paths are replaced by `${REPO}`, `${DATA}` and `${HOME}`.
