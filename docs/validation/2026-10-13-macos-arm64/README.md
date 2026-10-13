# ARM64 chaining and HG002 ONT WGS validation

The two ARM64 chaining changes preserve the tested Apple Silicon outputs and reduce
observed full-run time in separate nominal-1× experiments. The final build also
completed a nominal-10× HG002 ONT run to sorted, indexed BAM in **40.66 minutes**.
Native Linux ARM64/x86 WGS equality and independent alignment accuracy were not
established by these experiments.

This report concerns bonefish at source base
`2061243d0010e2895453a170c988ee54b37e7930` plus the changes in this PR.
The older rammap measurements in `docs/performance.md` describe a different tool
and workload. All measurements below were made on 2026-10-13.

## Implementation

`bonefish-core/src/align/chain_simd.rs` adds two NEON strategies already used by
the AVX2 backend:

- On dense anchor sets, merge predecessor lists from neighboring diagonal buckets
  in descending order, then score gathered candidates with the existing kernel.
  Negative signed reference/query coordinates retain the contiguous fallback;
  the bucket-density calculation uses i64 to avoid signed-coordinate overflow.
- On ordinary reads, reject off-band four-lane groups before loading the other
  fields, skip arithmetic for groups with no valid lane, and visit valid mask
  bits in descending order. The first chunk has 32 candidates, followed by chunks
  of 64. Masked and contiguous scans share const-specialized score arithmetic.
  Wrapping i32 diagonals retain the legacy signed-boundary behavior.

Both preserve predecessor order, tie selection, skip/visited bookkeeping,
no-skip behavior, scoring, backtracking and scratch reuse in the differential
fixtures. No preset or DP implementation changed. Production x86 kernels are
unchanged; the shared masked-score fixtures run against AVX2 on x86 and NEON on
ARM64. ARM-specific fixtures compare the optimized and legacy scan, including
DP scores and predecessor choices; scalar comparison is used where skip
semantics agree.

## System and input identity

| Property | Value |
|---|---|
| Host | Apple M4 Pro, 8 performance + 4 efficiency CPUs, 24 GiB RAM |
| Rust | rustc 1.98.1, aarch64-apple-darwin, LLVM 22.1.8 |
| Build | Generic ARM64, jemalloc, opt-level 3, fat LTO, one codegen unit; no native CPU flags |
| Other tools | minimap2 2.30-r1287; samtools 1.23.1 |
| Reads | Public HG002 ONT R10.4.1 ultra-long; Dorado 0.4.0+0aaf16d |
| Model | `dna_r10.4.1_e8.2_400bps_sup@v4.1.0` |
| Reference | Full GRCh38 no-alt analysis set: 195 contigs, 3,099,922,541 bases |
| Nominal-depth denominator | chr1–22, X and Y: 3,088,269,832 bases, including Ns |

Both datasets are **equal read-base prefixes of four source runs**, retaining
complete primary records in source order without length/quality filtering.
They are not uniform random samples. FASTQ retains sequence and quality but
discards the source uBAM modification tags. HG002 identity is from public source
provenance; these runs do not independently qualify donor identity.

The [source configuration](../../../scripts/validation/hg002_r1041_sources.json)
records all four public BAM URLs, lengths, multipart ETags and modification dates.
The 10× preflight rechecked those headers against the 1× sources. Multipart
ETags are source identity checks, not full-file SHA256 checksums. The complete
source BAMs were streamed only until each prefix target was reached.

| Input | Reads | Sequenced bases | Nominal depth | Read N50 |
|---|---:|---:|---:|---:|
| 1× | 86,553 | 3,088,602,282 | 1.000108× | 103,672 bp |
| 10× | 825,955 | 30,883,069,835 | 10.000120× | 108,853 bp |

The compressed reference SHA256 is
`3c8def6d325c5d1e934b2dd530c4d1709f027677a15859467071fddf3fff2026`.
The original combined 1× FASTQ SHA256 is
`960490bd561a62a82551a729abd97ed5f8fb27a228287b44c9fb6378b65b8116`.
Individual prefix hashes, counts, source headers and sampling receipts are in
[1× provenance](evidence/hg002-1x/DATASET_PROVENANCE.json) and
[10× provenance](evidence/hg002-10x/DATASET_PROVENANCE.json).

## Measured performance

| Comparison | Threads | Control | Candidate | Trials per build | Scope |
|---|---:|---:|---:|---:|---|
| Diagonal buckets | 8 | 142.198 s | 124.333 s | 2 | Fresh index + mapping + SAM, internal SSD |
| Masked scan, beyond buckets | 8 | 116.646 s | 106.888 s | 3 | Fresh index + mapping + SAM, internal SSD |
| Original release vs minimap2 | 4 | minimap2: 667.657 s | bonefish: 208.120 s | 1 | Fresh index + mapping + SAM; sorting excluded |
| Original release, native sorted 1× | 8 | — | 137.311 s | 1 | Index + mapping + sort + BAM + BAI, internal SSD |
| Final build, native sorted 10× | 8 | — | 2,439.845 s | 1 | Index + mapping + sort + BAM + BAI, external USB exFAT SSD |

The guarded bucket comparison reduced mean elapsed time by **12.6%** and summed
worker chaining duration by **39.2%**. Its final trials were `control-2`,
`candidate-3`, `candidate-4`, `control-3`, in that order. Earlier candidate-1/2
receipts describe the exploratory port before its signed-coordinate guard;
they and the earlier control-1 pilot are retained but excluded from final
averages. A differential test caught the boundary problem before the final
build. See the [summary](evidence/diagonal-buckets/BENCHMARK_RESULTS.json),
[metadata](evidence/diagonal-buckets/EXPERIMENT_METADATA.json) and individual logs.

The masked-scan comparison used the successful bucket build as control. Its
order was `control-1`, `candidate-1`, `candidate-2`, `control-2`, `control-3`,
`candidate-3`. Mean elapsed time fell a further **8.4%**, summed chaining
duration **30.2%**, and retired instructions **9.7%**. Control elapsed range
was 110.867–121.984 s; candidate range was 101.805–112.393 s. The ranges overlap.
These are small, separate trial cohorts on one machine; do not sum the
percentage gains or treat them as a general speed guarantee. See the
[summary](evidence/masked-scan/BENCHMARK_RESULTS.json).

minimap2 was measured on the smoke and 1× datasets, **not on 10×**. The 4-thread
1× wall-time ratio was 3.21×, using the original bonefish release and minimap2
`-ax lr:hq`; the outputs differ. No extrapolated 10× minimap2 time is presented
as a measurement. The optimized 8-thread 1× SAM trials and the 10× sorted-BAM
run have different output/storage scopes and cannot be compared as an aligner
speed ratio.

## Output validation

The original macOS release passed the repository's synthetic SAM/PAF variants,
aliases, raw-index equivalence, BAM/FASTQ input equivalence, sorted-BAM record
comparison and BAI query checks. Those checks are recorded in
[release-smoke](evidence/release-smoke/MACOS_TEST_RESULTS.json).

All four final bucket trials and all six masked-scan trials matched the
original 1× SAM exactly after removing `@PG`, with normalized SHA256
`4004e1b32a5515b8d006fa16dec0323f56af1c4938cabfbfbb3f0805829d4c66`.
Each had 132,946 records, 86,553 unique primary reads, 77,376 mapped primaries,
zero missing/extra primary names, and no warnings.

The 10× output passed complete BAM decoding, coordinate order, reference header
and bounds, CIGAR lengths, exact input primary inventory, and orientation-restored
sequence/quality BLAKE2b-128 fingerprints. Its BAI query passed. Every one of
the original 86,553 1× reads retained its exact decoded SAM record multiset:
**132,946 records compared**, including secondary/supplementary records,
headers excluded. See [BAM validation](evidence/hg002-10x/BAM_VALIDATION.json).

| 10× output | Count |
|---|---:|
| Total records | 1,285,384 |
| Primary reads | 825,955 |
| Mapped primary reads | 734,325 (88.9062%) |
| Unmapped primary reads | 91,630 |
| Secondary records | 218,649 |
| Supplementary records | 240,780 |

Primary-only nuclear mean depth was **9.0196×**, covering **94.1549%** of the
nuclear denominator. `samtools coverage --ff 0xF04 -q 0 -Q 0` includes MAPQ 0
and uses full contig lengths including Ns. Nominal sequenced depth and measured
mapped depth are different quantities. Per-contig values are retained in
[coverage](evidence/hg002-10x/PRIMARY_COVERAGE.tsv).

The 10× BAM SHA256 is
`a2cfa7503c3763ac60690a0a5e15a715779c9ade3b6785ce361e0880fdebf12c`;
its size is 27,450,389,542 bytes. The BAM, BAI and four FASTQs remain outside
Git; their locations/hashes are in [results](evidence/hg002-10x/RESULTS.json).

Fresh PR validation passed **129 core tests, 5 integration tests and 6 doctests**
(140 total); three existing benchmark tests remain ignored. Shared score
fixtures cover 1,431 chunks of lengths 0–64, unaligned starts, mixed targets and
signed-coordinate crossings. New Python validation checks cover corrupt inputs,
read inventory, reverse-strand qualities, census framing and RSS units. Small
real-tool smoke checks exercised the portable SAM runner, complete native-BAM
validator, subset comparison, and whole-record acquisition from all four sources.
Receipts are in [PR checks](evidence/pr-checks/VALIDATION.json).

## Remaining speed work and qualification limits

After the chaining changes, alignment occupied 48.45% and seeding 21.59% of
summed worker-stage wall durations in the three final 1× trials. These are not
serial elapsed times or exact hardware CPU shares.

The existing `census` feature recorded 7,241,438 helper calls and 36,208 sampled
input frames in one full 1× run. Recorded helper wall time totaled 291.670 s:
gap fill 138.508 s, left/right extensions together 133.847 s, gap reruns
19.315 s. The 117,812 gap helpers with maximum input length ≥1,024 consumed
75.527 s; the 5,923,275 helpers of length 128–255 consumed 48.147 s. Call count
alone does not identify the biggest time target. The portable decoder was also
rerun against the retained full traces and reproduced the saved summary.

Census timing includes helper probes/retries/subanchors, omits one equal-length
ungapped-first branch, and records outer bandwidth rather than each inner
kernel's bandwidth. Instrumented elapsed time is not a speed result. Sampled
helper inputs need reconstructed helper/probe arguments for faithful replay;
global sample IDs are not unique under concurrency. See
[census summary](evidence/dp-census/CENSUS_SUMMARY.json) and
[audit](evidence/dp-census/RESULT.txt). The approximately 505 MiB raw traces stay
outside Git, with sizes and hashes in the summary.

The default preset requests clean/adaptive DP, but native x86 AVX2 and legacy
ARM NEON already use different band semantics: scalar-exact bands on the clean
AVX2 path versus rounded vector bands/stale scores on legacy NEON; ARM adaptive
dispatch uses the static fallback. Porting the clean/adaptive path can change
scores, CIGARs and downstream choices. No DP change is included here. Real-data
x86/ARM equality must be established before enabling further such changes.

For 10×, the final sort merge/BAM/index phase took **1,073.362 s (17.89 min)**;
the preceding index/mapping/overlapping-spill phase took 1,366.213 s. There were
105 spill runs. Peak RSS was **13.974 GiB**; BSD time's separate peak-memory-footprint
measure was **21.811 GiB**, with zero swaps reported. A one-second stack sample
(included in the timer) saw 105 run-reader plus 105 decode threads, mostly
waiting. Device-wide I/O observations do not establish a process-specific root
cause. Bounded merge concurrency and controlled storage comparisons are future
experiments; output already has a 4-MiB buffer. No sorter change is included.

Preservation is not accuracy qualification. The 1× comparison retained **395
primary chromosome/strand disagreements** with minimap2, plus the smoke review
of a high-MAPQ placement disagreement. No new GIAB truth/variant-calling
evaluation was performed. Neither the Mac runs nor green small-fixture CI prove
cross-platform equality for this WGS workload. Native CI covers x86 Linux,
ARM Linux and Apple Silicon; real WGS parity and native performance remain
separate follow-up qualifications.

## Reproduction and artifact integrity

See [validation scripts](../../../scripts/validation/README.md) for build,
paired SAM trials, public-prefix acquisition, complete sorted-BAM validation and
census commands. Baseline/bucket/final build hashes and historical source patch
snapshots are retained with the experiment metadata. The final production
candidate SHA256 was
`1cb876c7cef4d1e3bebd433f9a6bfeea9e145669638687e89107a242004e669b`.

The evidence directory imports the small receipts, logs, placement differences,
source patch snapshots and audit reports. Old `committed: false`, `pushed: false`,
git-status strings and “uncommitted” report text describe the measurement time,
before this PR. They are historical evidence, not current branch status.

Only local paths were replaced: `${SESSION}` is the original session directory,
`${REPO}` its isolated checkout, `${REFERENCE}` the verified FASTA,
`${DATA_10X}` the external dataset directory, `${TOOLS}` the tool directory,
`${USER_HOME}` the local home path and `${EXTERNAL_SSD}` the external drive.
Public provider headers are preserved. The
[import manifest](evidence/IMPORT_MANIFEST.json) lists original and packaged
SHA256 values for each imported artifact; [bundle checksums](evidence/SHA256SUMS.txt)
cover the packaged evidence. The consistency checker verifies hashes and
recomputes trial means/reductions; it does not reread the external WGS BAM/FASTQs
or independently rerun mapping. Unrelated local dataset inventory is excluded.
