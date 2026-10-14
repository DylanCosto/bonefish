# Clean NEON/SSE4.1 DP: native validation and HG002 ONT WGS

On the Apple M4 Pro, the paired nominal-1× experiment reduced mean elapsed time
from **146.492 s to 121.805 s (16.85%)**, with three trials per build. Native
DP fixtures passed on x86 Linux, ARM Linux and Apple Silicon. The new clean DP
changes older ARM alignments: **3,538 reads had different SAM record sets**,
including **124 primary chromosome/strand placement changes**. Full real-WGS
x86/ARM equality and independent alignment accuracy were not tested here.

This report records [PR #2](https://github.com/DylanCosto/bonefish/pull/2), branch
`neon-clean-dp`, at exact tested source head
[`c1b7eb074b49806a86bbd0e48c617957ac436618`](https://github.com/DylanCosto/bonefish/commit/c1b7eb074b49806a86bbd0e48c617957ac436618),
against main
[`e8034d93ad7d7d213979147d648baaf98623117c`](https://github.com/DylanCosto/bonefish/commit/e8034d93ad7d7d213979147d648baaf98623117c).
Testing did not change aligner source. This later results commit publishes the
evidence; its commit ID is separate from the measured implementation.

The subsequent [HG001 replay](../2026-10-13-hg001-macos-parity/README.md) used the
user-supplied genome-wide sample and reproduced the supplied x86 MD5 prefix after
excluding only `@PG`. It records the full Mac hashes; the complete expected x86
hash and original command/reference receipt remain unprovided.

## Native and local checks

[Native CI run 38478508116](https://github.com/DylanCosto/bonefish/actions/runs/38478508116)
passed all three jobs at the exact tested head:

- x86 Linux, `ubuntu-24.04`: SSE4.1/AVX2 differential comparisons passed.
- ARM Linux, `ubuntu-24.04-arm`: clean NEON fixtures passed.
- Apple Silicon, `macos-14`: clean NEON fixtures passed.

The common 3,000-case result/CIGAR digest was **`0x4fcd48c4d9270cf3`**.
Cases cover deterministic random read-like inputs with Ns and indels, static,
adaptive and shifted bands, widths around 16/32 lanes, and scalar static-band
comparisons. Native fixture agreement is separate from a complete WGS replay.

Local validation passed **142 Rust tests** (131 core, five integration, six
doctests), with three existing benchmark tests ignored, and **13 Python checks**.
Existing evidence consistency checks also passed. End-to-end fixtures cover
SAM/PAF options, raw-index equivalence, BAM/FASTQ/BGZF input, and sorted BAM/BAI.
See [Rust log](evidence/cargo-tests.log), [Python log](evidence/python-tests.log),
[build log](evidence/build.log), [native CI log](evidence/native-ci.log) and
[machine-readable results](evidence/BENCHMARK_RESULTS.json).

## Input, build and timing identity

| Property | Value |
|---|---|
| Host | Apple M4 Pro, 12 CPUs (8 performance + 4 efficiency), 24 GiB RAM |
| Compiler | rustc 1.98.1, aarch64-apple-darwin, LLVM 22.1.8 |
| Build | Generic ARM64 release, jemalloc, opt-level 3, fat LTO, one codegen unit; no native CPU flags |
| Threads | 8 |
| Input | HG002 ONT R10.4.1 ultra-long, Dorado 0.4.0+0aaf16d SUP 4.1.0 |
| Reads / bases | 86,553 / 3,088,602,282; read N50 103,672 bp |
| Reference | GRCh38 no-alt analysis set, 195 contigs / 3,099,922,541 bases |
| Storage | FASTQ and SAM: internal APFS SSD; reference: external SSD for both builds |
| Timed work | Fresh FASTA-built index, mapping and SAM writing |
| Outside timing | Acquisition, input hashing/inventory, sequence/quality fingerprints, all output checks |

The dataset consists of equal read-base, whole-primary-record prefixes of four
public source runs in file order. It is not a uniform random sample. FASTQ retains
sequence and quality and discards source modification tags. The public source
URLs, headers, extraction method and part receipts are recorded in the existing
[1× provenance](../2026-10-13-macos-arm64/evidence/hg002-1x/DATASET_PROVENANCE.json)
and [source configuration](../../../scripts/validation/hg002_r1041_sources.json).
The reference and FASTQ hashes were rechecked for every trial.

| Artifact | SHA256 |
|---|---|
| Compressed FASTQ | `960490bd561a62a82551a729abd97ed5f8fb27a228287b44c9fb6378b65b8116` |
| Compressed reference | `3c8def6d325c5d1e934b2dd530c4d1709f027677a15859467071fddf3fff2026` |
| Main/control executable | `1cb876c7cef4d1e3bebd433f9a6bfeea9e145669638687e89107a242004e669b` |
| Branch/candidate executable | `40bfeb7053a0da22d1ce2a7709e77d441b9ca719afc684e949c38fd211109662` |

The control uses the previous optimized chaining implementation on main with
older ARM DP. The candidate adds clean NEON DP. Both executables use the same
compiler/profile/features, input files, output format and destination filesystem.
[Experiment metadata](evidence/EXPERIMENT_METADATA.json) records build commands,
paths as documented placeholders, thread count and predeclared trial order.

## All six measured trials

Order: control-1, candidate-1, candidate-2, control-2, control-3, candidate-3.
Every trial is included; all completed without warnings.

| Trial | Elapsed seconds | Receipt | Log |
|---|---:|---|---|
| control-1 | 170.339760 | [JSON](evidence/control-1.json) | [log](evidence/control-1.log) |
| candidate-1 | 130.570224 | [JSON](evidence/candidate-1.json) | [log](evidence/candidate-1.log) |
| candidate-2 | 115.939774 | [JSON](evidence/candidate-2.json) | [log](evidence/candidate-2.log) |
| control-2 | 132.658104 | [JSON](evidence/control-2.json) | [log](evidence/control-2.log) |
| control-3 | 136.477517 | [JSON](evidence/control-3.json) | [log](evidence/control-3.log) |
| candidate-3 | 118.903555 | [JSON](evidence/candidate-3.json) | [log](evidence/candidate-3.log) |

| Build | Mean elapsed | Range | Sample SD | Trials |
|---|---:|---:|---:|---:|
| Main / older ARM DP | 146.492 s | 132.658–170.340 s | 20.741 s | 3 |
| neon-clean-dp | 121.805 s | 115.940–130.570 s | 7.735 s | 3 |

Observed elapsed reduction is **16.85%**, a **1.203× throughput ratio**.
Mean summed worker alignment duration fell from 400.651 s to 299.671 s
(**25.20%**) and retired instructions fell **22.40%**. Worker durations sum
per-read wall times across concurrent work; they are not exact CPU shares or
serial elapsed time. Peak RSS and memory-footprint measurements are retained
per trial in [RUNS.json](evidence/RUNS.json).

These are three trials per build on one Mac with current memory, cache, thermal
and storage conditions. The first control trial is visibly slower and remains
included. No general speed guarantee, comparison to the earlier day's 106.9 s
mean, or 10×/30× timing claim follows from this experiment. This test writes SAM;
the earlier 10× experiment also sorted/compressed/indexed BAM on external storage
and has a different timing scope.

## Complete output checks and differences

All six runs passed exact primary read inventory and orientation-restored
sequence/quality BLAKE2b-128 checks. Each build produced the same normalized SAM
SHA256 in all three trials. The complete first control/candidate pair additionally
passed CIGAR full-query consumption, SEQ/CIGAR length, positive CIGAR operation,
reference bounds and identical reference-header checks for every record.

| Output property | Main/control | Branch/candidate |
|---|---:|---:|
| Records | 132,946 | 133,778 |
| Unique primary reads | 86,553 | 86,553 |
| Mapped primary reads | 77,376 | 77,376 |
| Unmapped primary reads | 9,177 | 9,177 |
| Missing / extra input primary names | 0 / 0 | 0 / 0 |

There were no lost or gained primary mappings. Record multisets were identical
for **83,015 reads** and changed for **3,538 reads**; **3,094 primary records**
changed. Of the changed fields, **124** both-mapped primary chromosome/strand
placements, **1,095** positions, **1,920** CIGARs and **682** MAPQs changed.
These categories overlap. Secondary/supplementary changes account for the net
832 additional records. [OUTPUT_COMPARISON.json](evidence/OUTPUT_COMPARISON.json)
contains all counts and examples, including 20 placement-change examples.
The [compressed changed-read JSONL](evidence/CHANGED_READS.jsonl.gz) retains
all 3,538 changed read IDs, record counts and before/after primary fields/CIGARs;
it does not contain every full SAM record.

Normalized SAM SHA256 excludes `@PG` and retains other headers and records:

- Control: `4004e1b32a5515b8d006fa16dec0323f56af1c4938cabfbfbb3f0805829d4c66`
- Candidate: `48da0e7d35654db9f6d3d45c171def634802c6ae4830b25c1fe06af54ff288d8`

The branch replaces legacy ARM DP with clean/adaptive semantics intended to
agree with x86 AVX2. **It does not preserve every older ARM alignment.**
The native fixtures agree across NEON/SSE4.1/AVX2, but a full WGS x86/ARM replay
and independent GIAB/variant-calling accuracy assessment remain unperformed.
The recorded candidate hash supplies an expected result for a native replay
using the exact same reference, reads and options.

Clean dispatch covers supported non-generic dual-affine scoring on NEON and
x86 AVX2/SSE4.1. CPUs lacking both x86 features, and generic-scoring paths,
retain legacy fallbacks; the report does not establish parity for those paths.

## Reproduction and packaged evidence

Build the recorded control and candidate commits in separate checkouts with the
same compiler and profile:

```bash
CARGO_PROFILE_RELEASE_STRIP=true cargo build --release --locked --features jemalloc
```

Use the existing [portable SAM runner](../../../scripts/validation/benchmark_sam.py)
for each of the six trials in the recorded order. Replace the binary, input and
scratch paths, and variant/trial values, for each run:

```bash
python3 -B scripts/validation/benchmark_sam.py /path/to/bonefish \
  --reference /path/to/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta.gz \
  --reads /path/to/HG002_R1041_UL.dorado0.4.0_sup4.1.0.nominal1x.fastq.gz \
  --threads 8 --variant candidate --trial 1 --output-dir /scratch/neon-clean-dp
```

For control runs, add `--expected-sam-sha256` with the control hash above.
Keep input/output on the same filesystems for both builds. The
[validation guide](../../../scripts/validation/README.md) describes acquisition,
input framing and failure handling. `BONEFISH_FORCE_SSE=1` selects SSE4.1 on
AVX2-capable x86 systems for comparison; native CI tested this dispatch.

Verify the packaged evidence without the large raw datasets:

```bash
python3 -B scripts/validation/verify_clean_dp_evidence.py
```

This hashes the bundle, recomputes all reported trial aggregates, checks receipt
and CI identity, and counts the compressed changed-read evidence. It does not
rerun WGS or establish biological accuracy. Local paths are replaced with
documented placeholders, and JSONL is losslessly compressed with deterministic
gzip. [IMPORT_MANIFEST.json](evidence/IMPORT_MANIFEST.json) records original and
packaged SHA256, sizes, transformations and placeholder meanings;
[SHA256SUMS.txt](evidence/SHA256SUMS.txt) covers every packaged evidence file.

The [analysis scripts](evidence/analysis/) are historical snapshots of the exact
local comparison, trial orchestration and summary calculation. Their path
placeholders must be replaced to run them; the comparator requires the complete,
validated HG002 SAMs in matching read-group order. `PR_RECEIPT.json`, `STATUS.json`
and `ORIGINAL_REPORT.md` retain measurement-time state before publication.
The executable, FASTQ and full SAMs remain outside Git. The first SAM pair is
retained in internal scratch; later SAMs were discarded only after complete
validation and recording their hashes.
