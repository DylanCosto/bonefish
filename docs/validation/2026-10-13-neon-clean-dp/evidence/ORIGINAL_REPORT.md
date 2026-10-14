# neon-clean-dp Mac validation

Branch: `neon-clean-dp`, head `c1b7eb074b49806a86bbd0e48c617957ac436618`.
Base: `e8034d93ad7d7d213979147d648baaf98623117c`; [draft PR #2](https://github.com/DylanCosto/bonefish/pull/2).
No branch source changes, commits or merges were made during testing.

Native CI passed on x86 Linux, ARM Linux and Apple Silicon. The SSE4.1/AVX2
comparisons and shared 3,000-case result/CIGAR digest passed. Local Mac validation
passed 142 Rust tests (131 core, 5 integration, 6 doctests), with 3 existing ignored
benchmarks, plus 13 Python checks and the previous receipt-consistency check.

The paired nominal-1x test used 86,553 HG002 ONT reads / 3,088,602,282 bases,
Dorado 0.4.0 SUP 4.1.0, equal-base file-order prefixes of four public source parts
against the complete GRCh38 no-alt analysis set. The sampling is not uniform random.
Both builds use rustc 1.98.1, generic ARM64 release settings, jemalloc, opt-level 3,
fat LTO, one codegen unit and 8 threads. Reference is on the external SSD; FASTQ
and SAM output are on the internal APFS SSD. Both use the same input/output paths
and command flags. Acquisition, input fingerprints/hashes and output validation
are outside the timer. Full recorded trial order: [['control', 1], ['candidate', 1], ['candidate', 2], ['control', 2], ['control', 3], ['candidate', 3]].

| Build | Mean elapsed | Range | Trials |
|---|---:|---:|---:|
| Current main / old ARM DP | 146.492 s | 132.658–170.340 s | 3 |
| neon-clean-dp | 121.805 s | 115.940–130.570 s | 3 |

Observed elapsed reduction: **16.85%**;
throughput ratio **1.203x**.
Summed worker alignment duration reduction: **25.20%**;
retired instruction reduction: **22.40%**.
These are current paired results, not a comparison against the earlier day's
106.9-second mean. Variation and memory/storage conditions limit generalization.

All six runs accounted for every input primary and preserved its sequence and
quality after restoring orientation. Outputs were deterministic within each build.
Control normalized SAM SHA256: `4004e1b32a5515b8d006fa16dec0323f56af1c4938cabfbfbb3f0805829d4c66`.
Candidate normalized SAM SHA256: `48da0e7d35654db9f6d3d45c171def634802c6ae4830b25c1fe06af54ff288d8`.

The complete first pair was compared as SHA256 multisets of decoded SAM records
per read, excluding headers. **3,538 reads**
have changed record sets; **3,094** changed primary
records. Both-mapped primary chromosome/strand changes: **124**.
Lost mapping: **0**; gained mapping: **0**.
Primary CIGAR changes: **1,920**;
MAPQ changes: **682**.

The branch replaces legacy ARM DP with clean/adaptive semantics intended to agree
with x86 AVX2. It does not preserve every older ARM alignment. Matching the native
kernel fixtures is separate from a complete real-WGS x86/ARM parity replay or an
independent accuracy assessment. Neither was performed here.

Clean dispatch covers NEON and AVX2/SSE4.1 for supported non-generic dual-affine
scoring. CPUs lacking both x86 features, and generic-scoring paths, retain the
legacy fallback. This qualification does not cover every custom configuration.

Detailed receipts, all trial logs, full output comparison and changed-read summaries
are retained here. The first control/candidate SAM pair is retained in internal
scratch; subsequent SAMs were discarded only after full checks and recording their
hashes, to bound disk use. Source stays at the exact tested head.
