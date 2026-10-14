#!/usr/bin/env python3
"""Verify archived clean-DP receipts, aggregates and changed-read evidence.

This checks recorded evidence, not external FASTQ/SAM data or alignment accuracy.
Run without Python's -O option: assertions enforce receipt consistency.
"""
import argparse
import collections
import gzip
import hashlib
import json
import math
from pathlib import Path
import statistics

REPO = Path(__file__).resolve().parents[2]
TESTED_HEAD = 'c1b7eb074b49806a86bbd0e48c617957ac436618'
BASE = 'e8034d93ad7d7d213979147d648baaf98623117c'


def load(path):
    return json.loads(path.read_text())


def close(actual, expected):
    assert math.isclose(actual, expected, rel_tol=1e-12, abs_tol=1e-9), (actual, expected)


def verify(evidence):
    manifest = load(evidence / 'IMPORT_MANIFEST.json')
    assert manifest['tested_source_head'] == TESTED_HEAD and manifest['base_commit'] == BASE
    imported = set()
    for item in manifest['artifacts']:
        path = evidence / item['path']
        assert path.resolve().is_relative_to(evidence.resolve())
        assert item['path'] not in imported, 'Duplicate import entry'
        imported.add(item['path'])
        data = path.read_bytes()
        assert len(data) == item['bytes']
        assert hashlib.sha256(data).hexdigest() == item['packaged_sha256'], path
        if item.get('compression'):
            decoded = gzip.decompress(data)
            assert len(decoded) == item['uncompressed_bytes']
            assert hashlib.sha256(decoded).hexdigest() == item['uncompressed_sha256']
            assert len(decoded.splitlines()) == item['records']
            if not item['local_paths_redacted']:
                assert item['uncompressed_sha256'] == item['original_sha256']
    checksummed = set()
    for line in (evidence / 'SHA256SUMS.txt').read_text().splitlines():
        expected, relative = line.split('  ', 1)
        path = evidence / relative
        assert path.resolve().is_relative_to(evidence.resolve())
        assert hashlib.sha256(path.read_bytes()).hexdigest() == expected, path
        assert relative not in checksummed, 'Duplicate checksum entry'
        checksummed.add(relative)
    files = {p.relative_to(evidence).as_posix() for p in evidence.rglob('*') if p.is_file()}
    assert checksummed == files - {'SHA256SUMS.txt'}, 'Missing/extra bundle checksum entry'
    assert imported == files - {'SHA256SUMS.txt', 'IMPORT_MANIFEST.json'}

    metadata = load(evidence / 'EXPERIMENT_METADATA.json')
    results = load(evidence / 'BENCHMARK_RESULTS.json')
    runs = load(evidence / 'RUNS.json')
    comparison = load(evidence / 'OUTPUT_COMPARISON.json')
    assert metadata['head_commit'] == results['source_head'] == TESTED_HEAD
    assert metadata['base_commit'] == results['base_commit'] == BASE
    assert results['source_unchanged_during_testing']
    assert runs == results['runs']
    assert [r['run'] for r in runs] == [f'{variant}-{trial}' for variant, trial in metadata['order']]
    expected_inputs = {
        '3c8def6d325c5d1e934b2dd530c4d1709f027677a15859467071fddf3fff2026',
        '960490bd561a62a82551a729abd97ed5f8fb27a228287b44c9fb6378b65b8116',
    }
    means = {}
    for variant in ('control', 'candidate'):
        group = results[variant]
        trials = sorted((r for r in runs if r['run'].startswith(variant + '-')), key=lambda r: r['run'])
        assert len(trials) == group['n'] == 3
        assert [r['run'] for r in trials] == group['runs']
        for run in trials:
            assert run == load(evidence / (run['run'] + '.json'))
            assert run['binary_sha256'] == metadata[variant + '_sha256']
            assert set(run['input_sha256'].values()) == expected_inputs
            assert run['returncode'] == 0 and not run['warnings'] and run['all_checks_passed']
            assert run['threads'] == metadata['threads'] == 8
            assert not run['census_instrumented']
            assert run['sam_sha256_without_pg'] == group['normalized_sam_sha256']
            assert run['records'] == group['records']
            assert run['primary'] == run['unique_primary_names'] == group['primary'] == 86553
            assert run['mapped_primary'] == group['mapped_primary'] == 77376
            assert run['missing_primary_names'] == run['unexpected_primary_names'] == 0
            assert run['input_sequences_and_qualities_preserved']
            if variant == 'control':
                assert run['same_as_original_release_sam']
            else:
                # Candidate trials supplied no legacy expected hash to the runner.
                assert run['same_as_original_release_sam'] is None
                assert run['sam_sha256_without_pg'] != results['control']['normalized_sam_sha256']
        times = [r['elapsed_seconds_including_index_and_sam'] for r in trials]
        assert times == group['elapsed_seconds']
        means[variant] = statistics.mean(times)
        close(group['mean_elapsed_seconds'], means[variant])
        close(group['min_elapsed_seconds'], min(times))
        close(group['max_elapsed_seconds'], max(times))
        close(group['sample_sd_elapsed_seconds'], statistics.stdev(times))
        for stage, duration in group['mean_worker_stage_sums_seconds'].items():
            close(duration, statistics.mean(r['worker_stage_sums_seconds'][stage] for r in trials))
        close(group['mean_instructions_retired'], statistics.mean(r['instructions_retired'] for r in trials))
        for key in ('peak_rss_bytes', 'peak_memory_footprint_bytes'):
            assert group[key] == [r[key] for r in trials]
    close(results['elapsed_reduction_percent'], 100 * (1 - means['candidate'] / means['control']))
    close(results['throughput_ratio'], means['control'] / means['candidate'])
    close(results['worker_alignment_duration_reduction_percent'], 100 * (1 -
          results['candidate']['mean_worker_stage_sums_seconds']['Alignment (Ext)'] /
          results['control']['mean_worker_stage_sums_seconds']['Alignment (Ext)']))
    close(results['retired_instructions_reduction_percent'], 100 * (1 -
          results['candidate']['mean_instructions_retired'] / results['control']['mean_instructions_retired']))

    assert comparison == results['output_comparison']
    counts = comparison['counts']
    assert all(comparison['validation'].values())
    assert counts['reads_compared'] == counts['same_record_multiset_reads'] + counts['changed_record_multiset_reads'] == 86553
    assert counts['both_mapped'] == 77376 and counts['both_unmapped'] == 9177
    assert counts.get('lost_mapping', 0) == counts.get('gained_mapping', 0) == 0
    assert counts['control_records'] == results['control']['records'] == 132946
    assert counts['candidate_records'] == results['candidate']['records'] == 133778
    observed = collections.Counter()
    read_ids = set()
    with gzip.open(evidence / 'CHANGED_READS.jsonl.gz', 'rt') as source:
        for line in source:
            row = json.loads(line)
            assert row['read_id'] not in read_ids, 'Duplicate changed-read ID'
            read_ids.add(row['read_id'])
            a, b = row['control_primary'], row['candidate_primary']
            observed['record_count_changed_reads'] += row['control_records'] != row['candidate_records']
            observed['record_delta'] += row['candidate_records'] - row['control_records']
            for field in ('flag', 'contig', 'position', 'mapq', 'cigar'):
                observed['primary_' + field + '_changed'] += a[field] != b[field]
            strand_change = bool(a['flag'] & 16) != bool(b['flag'] & 16)
            observed['primary_strand_changed'] += strand_change
            observed['both_mapped_primary_contig_or_strand_changed'] += (
                not (a['flag'] & 4 or b['flag'] & 4) and (a['contig'] != b['contig'] or strand_change))
    assert len(read_ids) == counts['changed_record_multiset_reads'] == 3538
    assert observed.pop('record_delta') == counts['candidate_records'] - counts['control_records']
    for key, value in observed.items():
        assert counts[key] == value, (key, counts[key], value)

    ci = results['native_ci']
    assert ci['headSha'] == TESTED_HEAD and ci['conclusion'] == 'success' and ci['status'] == 'completed'
    assert {j['name'] for j in ci['jobs']} == {'test (ubuntu-24.04)', 'test (ubuntu-24.04-arm)', 'test (macos-14)'}
    assert all(j['conclusion'] == 'success' and j['status'] == 'completed' for j in ci['jobs'])
    assert results['local_tests']['total_rust_passed'] == 142 and results['local_tests']['python_passed'] == 13
    print(f'Verified {len(imported)} imported artifacts, all bundle hashes, six trial aggregates, '
          f'{len(read_ids):,} changed-read entries and native CI identity at {TESTED_HEAD[:12]}.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--evidence', type=Path, default=REPO / 'docs/validation/2026-10-13-neon-clean-dp/evidence')
    verify(parser.parse_args().evidence)
