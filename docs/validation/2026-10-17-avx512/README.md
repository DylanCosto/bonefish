# AVX-512 alignment tuning and parity

The retained implementation speeds up wide dual-affine DP, chaining, sorting,
and minimizer sketching while preserving the tested SAM output. All measurements
below were made on 2026-10-17 in an x86_64 Linux VMware guest exposing an Intel
Xeon Gold 5418Y. These are historical measurements of the source and executables
identified in `evidence/build.json`, before PR packaging and integration with
subsequent documentation/packaging changes on `main`.

## Retained implementation

- Runtime dispatch uses the 64-lane clean dual-affine kernel for eligible bands
  wider than 32 cells. Narrow bands retain the 256-bit path; additional register
  paths use AVX-512 VBMI/VL only after checking those CPU features.
- The wide kernel uses shifted gap-state loads, mask-plane traceback, vectorized
  score tracking, and guarded gap arithmetic. Other scoring/band combinations
  retain the original arithmetic. Complete result and CIGAR fields are checked.
- Chaining uses masked predecessor scoring and vector merging without changing
  predecessor visit order. Shared work also benefits from stable counted RMQ
  ranks, dense range bounds, internal integer hashing, and radix-sort changes.
- Eligible sorting leaves use an AVX-512 network that preserves payload and tie
  order. The minimizer sketcher runs eight 64-bit streams, caches encoded bases,
  and reuses blocked prefix/suffix minima. The original emission order and
  duplicate/tie behavior are preserved.
- Existing runtime CPU detection and SIMD caps remain in force. AVX2, SSE and
  non-x86 fallbacks remain available. This series includes shared scalar/RMQ/radix
  and earlier narrow AVX2 recurrence improvements; it does not claim every
  fallback implementation is byte-for-byte unchanged from v0.1.6.

Singletrack, explicit-unroll variants, and subsequent experimental chaining
changes were not selected. No global `target-cpu=native`, global AVX-512 target
flags, PGO, or bootstrap compiler flags were used for the measured release build.

## Full-input parity and latest mapper comparison

The final minimizer increment was compared with the previous retained build,
which already contained the selected DP, chaining, and sorting changes. It is
not a comparison against unmodified v0.1.6 or minimap2.

All runs used 24 threads, the same prebuilt raw GRCh38 no-alt index (`k=19`,
`w=19`), and the complete original dataset: 161,762 reads, 2,713,864,134 bases.
The reference has 195 contigs and 3,099,922,541 bases. Three alternating timing
pairs produced these medians:

| Metric | Previous retained build | Selected build |
| --- | ---: | ---: |
| Elapsed seconds | 32.86 | 31.84 |
| User + system CPU seconds | 550.12 | 536.36 |
| Peak RSS, GiB | 13.815 | 13.808 |

This is an observed 1.0320x elapsed throughput improvement and 2.501% less CPU
time. Each paired observation reduces CPU time. Other user setup jobs overlapped
the measurements, so elapsed-time gains remain provisional. Index construction,
input warming, external coordinate sorting, and BAM compression are excluded;
index loading, decoding, mapping, and SAM emission are included.

All 11 full-input checks passed: the six timed runs, forced AVX2 on both builds,
forced SSE on the selected build, the disabled-minimizer path, and a separate
activation run. Canonicalization removes only the command-dependent `@PG` header.
Every result has SAM MD5 `9bc7ff3e9f498cc1530b2b1d78c47e10`. Activation counters
recorded 143,152 new-kernel calls and 2,706,881,735 aggregate bases; those counters
were disabled for timing. Input SHA-256 values and all individual observations
are in `evidence/mapper.json`; original GNU-time, checksum, and mapper-log receipts
are included alongside it.

## Kernel measurements are separate from mapper throughput

The retained wide DP implementation previously measured **1.895913x versus
AVX2** on 805 captured calls (66 exact-tracking and 739 approximate-tracking),
including allocation, initialization, score tracking, and CIGAR traceback.
Six alternating outer batches each contained three alternating trials with
eight corpus repetitions on pinned CPU 23. An actual-checkout confirmation
measured 1.882631x. Those kernel sources are retained in the final build.
`evidence/wide-dp-report.txt`, `wide-dp-kernel-replay.csv`, and `wide-dp-proof.txt`
preserve that earlier stage's results and arithmetic conditions. References to
the selected executable or commit status in that report describe that stage.

The final blocked minimizer measured **1.7246x versus AVX2** in a separate
diagnostic replay of the first 256 original reads (4,571,254 bases), using the
actual index parameters and CPU 23. All minimizer fields and their order matched.
The replay reports thread CPU time, not physical cycles; its receipts are in
`evidence/minimizer-replay.json` and `minimizer-replay.csv`.

Neither kernel ratio is a whole-mapper ratio. These measurements do not establish
a 2x whole-mapper improvement or a cumulative 2.06x speedup. Earlier full 10x
DP-stage parity results are recorded in the historical wide-DP report; the final
minimizer increment was validated on the complete original 1x dataset.

## Tidy before merge

Before merging, the series was cut down to the code its default configuration runs: 32
diagnostic `BONEFISH_*` switches, seven non-default chaining variants, the DP shifted-load
and traceback alternatives, the replay benchmarks and the minimizer activation counter were
removed, and the AVX-512 sorting network now follows the process SIMD cap like the other
AVX-512 paths. AVX-512 tests compare against the AVX2 or scalar code; CI lists the runner's
SIMD features in `BONEFISH_REQUIRE_SIMD`, which turns a skipped test into a failure. The
measurements above were made before this step.

On an Intel i5-12600K (AVX2, no AVX-512), 16 threads, raw GRCh38 index, the tidied build
gives the same output as v0.1.6 (SAM MD5 without `@PG`):

| reads | v0.1.6 = tidied build |
|---|---|
| genome-wide HG001 sample (default, `BONEFISH_FORCE_SSE=1`, `BONEFISH_FORCE_SCALAR=1`) | `9bc7ff3e9f498cc1530b2b1d78c47e10` |
| HG002 chr22, PAW71238 (27x) | `7e381581bc6bb5a4d56ca655b4c278b1` |
| HG002 chr20, PAW70337 (55x) | `165b75dbfa1282faddd2c04cd86a4efb` |

Genome-wide sample, three alternating runs each: v0.1.6 672.7-677.0 CPU-s (44.7-47.1 s),
tidied build 650.8-653.1 CPU-s (45.7-48.0 s). The shared changes (RMQ ranks, radix sort,
chaining buckets) save about 3% of CPU time on AVX2 machines too.

On the same Xeon Gold 5418Y VM (24 threads), the tidied build (commit `ca4270d`) passed
`simd_parity.sh` on 20 Oct: default (AVX-512), forced AVX2 and forced SSE4.1 all give
`9bc7ff3e9f498cc1530b2b1d78c47e10`. Default used about 526 CPU-s (27.6-29.7 s wall), forced
AVX2 about 647 and SSE4.1 about 855. On that VM v0.1.6 used 721-723 CPU-s by default (its
AVX-512 path) and 671-685 with AVX2, so this series saves about 27% of CPU time on AVX-512
machines. A whole-program `-C target-cpu=x86-64-v4` build saved only a further 5 CPU-s
(520.3 against 525.1, three alternating rounds; wall time within noise), so release builds
stay portable with runtime dispatch.

## Tests and reproduction

The freshly built selected checkout passed **188 tests, zero failures, and seven
intentionally ignored external-data/benchmark tests**, with no diagnostic saying
AVX-512 was unavailable. The binary was byte-identical to the independently
tested private candidate. See `evidence/test-results.txt` and `build.json`.
The recorded source hashes identify all 16 selected alignment source files.
PR packaging removes an extra blank line at EOF in each register-kernel file;
`build.json` records both hashes. These edits change no code.

Build and run the normal release checks:

```bash
cargo test --release --locked --features jemalloc
cargo build --release --locked --features jemalloc --bin bonefish
python3 -B -m unittest discover -s scripts/validation -p 'test_*.py'
python3 -B scripts/validation/verify_evidence.py
python3 -B scripts/validation/verify_clean_dp_evidence.py
```

With the same reference and reads available, the repository's parity script
provides default/forced-AVX2/forced-SSE checks and alternating speed comparisons:

```bash
THREADS=24 bash scripts/validation/simd_parity.sh \
  target/release/bonefish \
  "$BENCH_ROOT/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz" \
  "$BENCH_ROOT/hg001_wg_1x.fq.gz" \
  "$BENCH_ROOT/avx512-pr-parity"
```

The measurements used diagnostic switches (`BONEFISH_FORCE_DP_AVX2`,
`BONEFISH_SKETCH_AVX512_OFF`, `BONEFISH_SKETCH_AVX512_REPORT` and others) that were
removed before merging (see "Tidy before merge"); the evidence files still name them. Use
the global SIMD force variables (`BONEFISH_FORCE_AVX2=1`, `BONEFISH_FORCE_SSE=1`) for
fallback checks.
Check the imported evidence files independently with:

```bash
cd docs/validation/2026-10-17-avx512
sha256sum -c SHA256SUMS
```

VMware exposed no usable cycle or APERF/MPERF counters, even to root, so actual
physical frequency and AVX-512 downclocking were not measured. Variant-calling
F1/Clair3 validation was not run; exact within-Bonefish SAM parity does not
establish cross-tool variant-calling accuracy. Cross-platform PR CI is separate
from the recorded VM results.
