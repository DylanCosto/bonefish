# Contributing

Bug reports, test data and pull requests are welcome. For bugs, please use the issue template:
the version, command line, CPU and a few reads that reproduce the problem make it fast to fix.

## Build and test

```bash
cargo build --release --features jemalloc
cargo test --release --features jemalloc       # unit tests, golden-output regression, doc tests
cargo clippy --release --features jemalloc --all-targets
```

CI runs the same on Linux x86_64, Linux arm64 and macOS (arm64). On x86 runners with AVX-512,
CI requires the AVX-512 tests to run rather than skip.

## Output changes

The golden-output test (`bonefish/tests/regression.rs`) pins the exact output of every preset on
small synthetic datasets in `bonefish/tests/data/`. Every change therefore falls into one of two
cases:

- **Speed or refactoring:** the output must stay byte-identical, and the regression test must pass
  unchanged. For larger changes, also compare a real dataset before and after (SAM with `@PG`
  lines removed). The validation reports quote MD5s for the genome-wide HG001 sample.
- **Intended output changes:** regenerate the expected hashes and explain the change in the pull
  request:

  ```bash
  BONEFISH_BLESS=1 cargo test --release --features jemalloc --test regression -- --nocapture
  ```

  Paste the printed rows into `bonefish/tests/data/expected.tsv`. An output change also needs
  evidence that accuracy does not get worse. For the `dna` preset, that means variant calling
  against GIAB (Clair3 for small variants, Sniffles2 for SVs) on at least one chromosome. For the
  spliced presets, it means junction accuracy against an annotation and against simulated reads
  with known intron chains. `docs/validation/` shows how earlier changes were validated, with
  scripts.

The test datasets are generated deterministically by `scripts/make_regression_data.py` (DNA) and
`scripts/make_regression_rna_data.py` (RNA/cDNA).

## Code

- **Style:** match the surrounding code.
- **SIMD:** kernels need a scalar or portable fallback, selected at run time. Tests that compare
  them with the reference implementation should skip through
  `crate::align::dp::simd_test_available` when the CPU lacks the feature.
- **DNA speed:** new presets or options must not slow the `dna` preset. Keep per-read and per-call
  branches out of the DP and chaining inner loops.
