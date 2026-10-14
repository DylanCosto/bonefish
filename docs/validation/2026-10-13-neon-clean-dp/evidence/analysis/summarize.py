import datetime
import hashlib
import json
from pathlib import Path
import statistics
import subprocess

ROOT = Path('${SESSION}')
OUT = ROOT / 'outputs/bonefish-neon-clean-dp'
REPO = ROOT / 'work/bonefish-neon-clean-dp/src'
config = json.loads((OUT / 'EXPERIMENT_METADATA.json').read_text())
runs = json.loads((OUT / 'RUNS.json').read_text())
comparison = json.loads((OUT / 'OUTPUT_COMPARISON.json').read_text())
assert [run['run'] for run in runs] == [f'{variant}-{trial}' for variant, trial in config['order']]
for run in runs:
    assert run['all_checks_passed'] and run['returncode'] == 0 and not run['warnings']
    assert run['primary'] == run['unique_primary_names'] == 86553
    assert run['missing_primary_names'] == run['unexpected_primary_names'] == 0
    assert run['input_sequences_and_qualities_preserved']
    assert run['input_sha256'][config['reads']] == '960490bd561a62a82551a729abd97ed5f8fb27a228287b44c9fb6378b65b8116'
    assert run['input_sha256'][config['reference']] == '3c8def6d325c5d1e934b2dd530c4d1709f027677a15859467071fddf3fff2026'

def aggregate(variant):
    rs = [run for run in runs if run['run'].startswith(variant)]
    times = [run['elapsed_seconds_including_index_and_sam'] for run in rs]
    assert len({run['sam_sha256_without_pg'] for run in rs}) == 1, 'Nondeterministic output'
    assert all(run['binary_sha256'] == config[variant + '_sha256'] for run in rs)
    stages = {key: statistics.mean(run['worker_stage_sums_seconds'][key] for run in rs)
              for key in rs[0]['worker_stage_sums_seconds']}
    return dict(n=len(rs), runs=[run['run'] for run in rs], elapsed_seconds=times,
                mean_elapsed_seconds=statistics.mean(times), min_elapsed_seconds=min(times), max_elapsed_seconds=max(times),
                sample_sd_elapsed_seconds=statistics.stdev(times), mean_worker_stage_sums_seconds=stages,
                peak_rss_bytes=[run['peak_rss_bytes'] for run in rs],
                peak_memory_footprint_bytes=[run['peak_memory_footprint_bytes'] for run in rs],
                mean_instructions_retired=statistics.mean(run['instructions_retired'] for run in rs),
                normalized_sam_sha256=rs[0]['sam_sha256_without_pg'],
                records=rs[0]['records'], primary=rs[0]['primary'], mapped_primary=rs[0]['mapped_primary'])

control, candidate = aggregate('control'), aggregate('candidate')
ci = json.loads(subprocess.check_output(['gh', 'run', 'view', '38478508116', '--repo', 'DylanCosto/bonefish',
                                         '--json', 'url,headSha,status,conclusion,jobs'], text=True))
assert ci['headSha'] == config['head_commit'] and ci['conclusion'] == 'success'
assert len(ci['jobs']) == 3 and all(job['conclusion'] == 'success' for job in ci['jobs'])
assert subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=REPO, text=True).strip() == config['head_commit']
assert not subprocess.check_output(['git', 'status', '--porcelain'], cwd=REPO, text=True).strip()
elapsed_reduction = 100 * (1 - candidate['mean_elapsed_seconds'] / control['mean_elapsed_seconds'])
ext_reduction = 100 * (1 - candidate['mean_worker_stage_sums_seconds']['Alignment (Ext)'] / control['mean_worker_stage_sums_seconds']['Alignment (Ext)'])
instruction_reduction = 100 * (1 - candidate['mean_instructions_retired'] / control['mean_instructions_retired'])
result = dict(control=control, candidate=candidate, elapsed_reduction_percent=elapsed_reduction,
              throughput_ratio=control['mean_elapsed_seconds'] / candidate['mean_elapsed_seconds'],
              worker_alignment_duration_reduction_percent=ext_reduction,
              retired_instructions_reduction_percent=instruction_reduction,
              output_comparison=comparison, runs=runs, native_ci=ci,
              source_head=config['head_commit'], base_commit=config['base_commit'],
              source_unchanged_during_testing=True, timing_scope=config['timing_scope'],
              hardware=dict(chip='Apple M4 Pro', cpus=12, ram_bytes=25769803776),
              local_tests=dict(core_passed=131, integration_passed=5, doctests_passed=6, total_rust_passed=142,
                               existing_ignored_benchmarks=3, python_passed=13),
              qualifications=['Three trials per build on one Mac; small sample and current system/storage conditions, no general performance guarantee.',
                              'Worker stage durations sum per-read wall times, not serial elapsed time or exact CPU shares.',
                              'Native clean-kernel parity is checked on 3,000 deterministic cases, not a full x86 WGS replay.',
                              'Old Mac alignments used legacy/static DP. Output changes to clean/adaptive behavior are reported separately from structural preservation and accuracy.',
                              'Clean dispatch covers NEON and x86 AVX2/SSE4.1 for supported non-generic dual-affine scoring. CPUs lacking both x86 features and generic-scoring paths retain legacy fallbacks.',
                              'No new independent GIAB/variant-calling accuracy qualification was performed.'],
              completed_utc=datetime.datetime.now(datetime.timezone.utc).isoformat())
(OUT / 'BENCHMARK_RESULTS.json').write_text(json.dumps(result, indent=2) + '\n')
counts = comparison['counts']
report = f'''# neon-clean-dp Mac validation

Branch: `neon-clean-dp`, head `{config['head_commit']}`.
Base: `{config['base_commit']}`; [draft PR #2](https://github.com/DylanCosto/bonefish/pull/2).
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
are outside the timer. Full recorded trial order: {config['order']}.

| Build | Mean elapsed | Range | Trials |
|---|---:|---:|---:|
| Current main / old ARM DP | {control['mean_elapsed_seconds']:.3f} s | {control['min_elapsed_seconds']:.3f}–{control['max_elapsed_seconds']:.3f} s | 3 |
| neon-clean-dp | {candidate['mean_elapsed_seconds']:.3f} s | {candidate['min_elapsed_seconds']:.3f}–{candidate['max_elapsed_seconds']:.3f} s | 3 |

Observed elapsed reduction: **{elapsed_reduction:.2f}%**;
throughput ratio **{result['throughput_ratio']:.3f}x**.
Summed worker alignment duration reduction: **{ext_reduction:.2f}%**;
retired instruction reduction: **{instruction_reduction:.2f}%**.
These are current paired results, not a comparison against the earlier day's
106.9-second mean. Variation and memory/storage conditions limit generalization.

All six runs accounted for every input primary and preserved its sequence and
quality after restoring orientation. Outputs were deterministic within each build.
Control normalized SAM SHA256: `{control['normalized_sam_sha256']}`.
Candidate normalized SAM SHA256: `{candidate['normalized_sam_sha256']}`.

The complete first pair was compared as SHA256 multisets of decoded SAM records
per read, excluding headers. **{counts['changed_record_multiset_reads']:,} reads**
have changed record sets; **{counts['primary_record_changed']:,}** changed primary
records. Both-mapped primary chromosome/strand changes: **{counts['both_mapped_primary_contig_or_strand_changed']:,}**.
Lost mapping: **{counts.get('lost_mapping',0):,}**; gained mapping: **{counts.get('gained_mapping',0):,}**.
Primary CIGAR changes: **{counts['primary_cigar_changed']:,}**;
MAPQ changes: **{counts['primary_mapq_changed']:,}**.

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
'''
(OUT / 'README.md').write_text(report)
files = sorted(p for p in OUT.iterdir() if p.is_file() and p.name not in ('SHA256SUMS.txt', 'STATUS.json'))
(OUT / 'SHA256SUMS.txt').write_text(''.join(f'{hashlib.sha256(path.read_bytes()).hexdigest()}  {path.name}\n' for path in files))
print(json.dumps({key: result[key] for key in ('elapsed_reduction_percent', 'worker_alignment_duration_reduction_percent', 'retired_instructions_reduction_percent')}, indent=2))
print(json.dumps(counts, indent=2))
