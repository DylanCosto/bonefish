# Reproduce the HG002 ONT validation

The [Mac report](../../docs/validation/2026-10-13-macos-arm64/README.md) explains
the measured builds, sampling, timing scope, output checks and limits. These
scripts are configurable versions of the local experiment runners; the report's
archived numbers are historical measurements, not timings of the Python fixture
checks. Requirements: Python 3.10+, Rust, samtools (tested 1.23.1), and BSD/GNU
`/usr/bin/time`. The source streamer needs samtools built with HTTPS support.

Run from the repository root. Keep generated reads/BAMs/traces and runner
directories outside Git. Do not run the validation scripts with `python -O`:
their assertions check malformed input and receipt consistency.

## Build and fast checks

```bash
CARGO_PROFILE_RELEASE_STRIP=true cargo build --release --locked --features jemalloc
CARGO_PROFILE_RELEASE_STRIP=true cargo test --release --locked --features jemalloc
python3 -B -m unittest discover -s scripts/validation -p 'test_*.py'
python3 -B scripts/validation/verify_evidence.py
```

The existing CI matrix runs the Rust tests, Python fixture checks and receipt
checker on x86 Linux, ARM Linux and Apple Silicon. It does not acquire/run human
WGS or prove native WGS output equality. The receipt checker hashes the imported
evidence and recomputes the final trial aggregates; it does not have access to
the external raw data.

For the original control, build base
`2061243d0010e2895453a170c988ee54b37e7930` in a separate checkout with the same
compiler/features/profile. For the intermediate bucket control, apply
`docs/validation/2026-10-13-macos-arm64/evidence/diagonal-buckets/neon-diagonal-buckets.patch`
to that base before building. The final source is this PR. Historical patches
are snapshots for reproduction; use the live source for development. Executable
SHA256 values are recorded in the experiment metadata; compiler changes may
change executable bytes even when behavior is unchanged.

## Acquire and validate a nominal 1×, 10× or 30× prefix

Supply the existing compressed GRCh38 no-alt analysis-set FASTA. Its required
SHA256 is recorded in `hg002_r1041_sources.json`; a different compressed reference
is rejected. Reference lengths are read from FASTA, so ordinary gzip works
without a samtools FASTA index. Acquisition checks the recorded public URLs,
sizes, ETags and dates before streaming whole primary records from the four
uBAMs. Prefix extraction discards modification tags; it is not random sampling.

```bash
python3 -B scripts/validation/run_hg002_wgs.py \
  --bonefish /path/to/bonefish \
  --reference /path/to/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta.gz \
  --coverage 10 --threads 8 --sort-mem 1G \
  --data-dir /scratch/hg002-10x/data \
  --work-dir /scratch/hg002-10x/work \
  --output-dir /scratch/hg002-10x/receipts
```

Change `--coverage` and directories to obtain 1× or 30×. The minimum-free-space
preflight defaults to 15 GiB per nominal x, with `--minimum-free-gib` available
for the actual storage plan; this estimate is not a disk-capacity guarantee.
The original 10× test used an external USB exFAT SSD for reads, spills and BAM.
Use the same output/storage configuration for any performance comparison.

To compare all records belonging to an earlier subset, add both options:

```bash
  --previous-bam /path/to/bonefish.1x.sorted.bam \
  --previous-reads /path/to/nominal1x.fastq.gz
```

The runner verifies gzip CRC/framing, unique read IDs, exact primary inventory,
orientation-restored sequence/quality fingerprints, complete BAM decode, CIGAR
and reference bounds, coordinate order, headers, BAI query and optional decoded
record multisets for the earlier subset. Primary nuclear coverage is measured
with `samtools coverage --ff 0xF04 -q 0 -Q 0`.

Index construction, mapping, overlapping sort spills, final merge/compression
and BAI generation are timed. Acquisition, input/output validation, coverage
and checksum calculation occur outside timing. BSD RSS is in bytes; GNU RSS
in KiB is converted to bytes. Unavailable footprint/instruction metrics remain
null. Failures keep generated evidence, incomplete FASTQs and partial BAMs.
Matching completed prefix receipts can be reused before alignment; existing
alignment/results are not overwritten. Never combine a receipt with changed
source data or an unrelated BAM.

## Paired SAM timing and output preservation

Use the same reference, FASTQs, compiler, generic target, features, threads and
output filesystem for control/candidate. Choose trial order beforehand and
retain every run, including failures/pilots. Warm-cache and thermal/background
conditions can change timing. The reported experiments' actual orders and
excluded exploratory builds are retained in their metadata.

```bash
python3 -B scripts/validation/benchmark_sam.py /path/to/control/bonefish \
  --reference /path/to/reference.fasta.gz --reads /path/to/nominal1x.fastq.gz \
  --threads 8 --variant control --trial 1 --output-dir /scratch/paired
python3 -B scripts/validation/benchmark_sam.py /path/to/candidate/bonefish \
  --reference /path/to/reference.fasta.gz --reads /path/to/nominal1x.fastq.gz \
  --threads 8 --variant candidate --trial 1 --output-dir /scratch/paired \
  --expected-sam-sha256 4004e1b32a5515b8d006fa16dec0323f56af1c4938cabfbfbb3f0805829d4c66
```

The expected hash above applies only to this exact recorded 1× input/reference
and SAM options. On a different dataset first record the control hash, then
pass it to subsequent trials. The hash includes headers except `@PG`; comparing
different pipeline formats requires decoded records instead. Each trial checks
all primary reads and sequence/quality preservation after timing. It retains
SAM, stderr/time log and JSON receipt; use new trial IDs to avoid overwriting.
Whole input and binary hashes are recorded outside timing. To reproduce the
original 1× FASTQ layout from newly acquired parts, concatenate the four gzip
members in part order into one `.fastq.gz`; compressed-file bytes may differ
from the archived artifact even when read records are identical.

## Census audit

Build with the existing feature and use a separate output/prefix:

```bash
CARGO_PROFILE_RELEASE_STRIP=true cargo build --release --locked --features jemalloc,census
python3 -B scripts/validation/benchmark_sam.py target/release/bonefish \
  --reference /path/to/reference.fasta.gz --reads /path/to/nominal1x.fastq.gz \
  --threads 8 --variant candidate --trial 1 --output-dir /scratch/census-run \
  --census-prefix /scratch/census/profile \
  --expected-sam-sha256 4004e1b32a5515b8d006fa16dec0323f56af1c4938cabfbfbb3f0805829d4c66
python3 -B scripts/validation/analyze_census.py /scratch/census/profile /scratch/census/summary.json
```

The runner explicitly sets `BONEFISH_CENSUS` and sampling every 200 calls, with
no minimum-length sampling override; this is the decoder's expected sampling
configuration. It refuses an existing census prefix. The decoder checks the
68-byte call records, 117-byte sample headers, sample/input lengths, nt4 alphabet,
sampling counts and CIGAR-consumption bounds. It summarizes helper wall durations,
not every inner kernel or exact CPU shares. Consult the report's replay/coverage
caveats before using samples for an optimization. Rebuild without `census`
before production timing; instrumented times are not speedup measurements.
