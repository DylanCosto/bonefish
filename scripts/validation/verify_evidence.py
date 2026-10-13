#!/usr/bin/env python3
"""Check packaged receipt integrity and recompute reported trial aggregates.

This checks the evidence bundle, not the external FASTQ/BAM data or mapping accuracy.
"""
import argparse
import hashlib
import json
import math
from pathlib import Path
import statistics

REPO = Path(__file__).resolve().parents[2]
SAM_HASH = '4004e1b32a5515b8d006fa16dec0323f56af1c4938cabfbfbb3f0805829d4c66'
CANDIDATE_HASH = '1cb876c7cef4d1e3bebd433f9a6bfeea9e145669638687e89107a242004e669b'


def load(path):
    return json.loads(path.read_text())


def close(actual, expected):
    assert math.isclose(actual, expected, rel_tol=1e-12, abs_tol=1e-9), (actual, expected)


def verify(evidence):
    manifest = load(evidence / 'IMPORT_MANIFEST.json')
    paths = set()
    for item in manifest['artifacts']:
        path = evidence / item['path']
        assert path.resolve().is_relative_to(evidence.resolve())
        assert item['path'] not in paths, 'Duplicate import entry'
        paths.add(item['path'])
        data = path.read_bytes()
        assert len(data) == item['bytes']
        assert hashlib.sha256(data).hexdigest() == item['packaged_sha256'], path
    checksummed = set()
    for line in (evidence / 'SHA256SUMS.txt').read_text().splitlines():
        expected, relative = line.split('  ', 1)
        path = evidence / relative
        assert path.resolve().is_relative_to(evidence.resolve())
        assert hashlib.sha256(path.read_bytes()).hexdigest() == expected, path
        assert relative not in checksummed, 'Duplicate checksum entry'
        checksummed.add(relative)
    assert checksummed == {path.relative_to(evidence).as_posix() for path in evidence.rglob('*')
                           if path.is_file() and path.name != 'SHA256SUMS.txt'}, 'Missing/extra bundle checksum entry'

    for directory, reduction_key in [('diagonal-buckets', 'elapsed_reduction_percent'),
                                      ('masked-scan', 'additional_elapsed_reduction_percent')]:
        folder = evidence / directory
        metadata = load(folder / 'EXPERIMENT_METADATA.json')
        summary = load(folder / 'BENCHMARK_RESULTS.json')
        means = {}
        for variant in ('control', 'candidate'):
            group = summary[variant]
            runs = [load(folder / f'{name}.json') for name in group['runs']]
            assert len(runs) == group['n']
            binary_hash = metadata['control_sha256'] if variant == 'control' else metadata.get('candidate_sha256', metadata.get('final_candidate_sha256'))
            for run in runs:
                assert run['binary_sha256'] == binary_hash
                assert run['returncode'] == 0 and not run['warnings']
                assert run['sam_sha256_without_pg'] == SAM_HASH
                assert run['same_as_original_release_sam']
                assert run['records'] == 132946
                assert run['primary'] == run['unique_primary_names'] == 86553
                assert run['mapped_primary'] == 77376
                assert run['missing_primary_names'] == run['unexpected_primary_names'] == 0
            times = [run['elapsed_seconds_including_index_and_sam'] for run in runs]
            assert times == group['elapsed_seconds']
            means[variant] = statistics.mean(times)
            close(group['mean_elapsed_seconds'], means[variant])
            close(group['min_elapsed_seconds'], min(times))
            close(group['max_elapsed_seconds'], max(times))
            close(group['mean_worker_chaining_seconds'], statistics.mean(run['worker_stage_sums_seconds']['Chaining'] for run in runs))
        close(summary[reduction_key], 100 * (1 - means['candidate'] / means['control']))
        close(summary['worker_chaining_duration_reduction_percent'], 100 * (1 - summary['candidate']['mean_worker_chaining_seconds'] / summary['control']['mean_worker_chaining_seconds']))

    original = load(evidence / 'hg002-1x/DATASET_PROVENANCE.json')
    large = load(evidence / 'hg002-10x/RESULTS.json')
    dataset, run, validation = large['dataset'], large['bonefish_run'], large['validation']
    assert original['reference_sha256'] == dataset['reference_sha256']
    assert dataset == load(evidence / 'hg002-10x/DATASET_PROVENANCE.json')
    assert run == load(evidence / 'hg002-10x/bonefish.run.json')
    assert validation == load(evidence / 'hg002-10x/BAM_VALIDATION.json')
    assert dataset['reads'] == 825955 and dataset['bases'] == 30883069835
    assert dataset['reads'] == sum(part['reads'] for part in dataset['source_parts'])
    assert dataset['bases'] == sum(part['bases'] for part in dataset['source_parts'])
    close(dataset['nominal_read_base_coverage'], dataset['bases'] / dataset['coverage_denominator_bases'])
    assert run['binary_sha256'] == CANDIDATE_HASH and run['exit_code'] == 0 and not run['warnings']
    assert validation['counts']['records'] == 1285384
    assert validation['counts']['primary'] == dataset['reads']
    assert validation['counts']['mapped_primary'] + validation['counts']['unmapped_primary'] == dataset['reads']
    assert sum(validation['counts'][key] for key in ('primary', 'secondary', 'supplementary')) == validation['counts']['records']
    assert all(validation[key] for key in ('quickcheck_passed', 'full_decode_passed', 'coordinate_sort_passed',
                                          'exact_primary_read_inventory', 'input_sequences_and_qualities_preserved',
                                          'cigar_and_reference_bounds_passed', 'reference_header_matches',
                                          'previous_1x_record_sets_identical'))
    assert validation['previous_1x_reads_compared'] == original['reads']
    assert validation['previous_1x_alignment_records_compared'] == 132946
    assert large['all_checks_passed']
    profile = load(evidence / 'dp-census/profile-1.json')
    census = load(evidence / 'dp-census/CENSUS_SUMMARY.json')
    assert profile['same_as_original_release_sam'] and profile['sam_sha256_without_pg'] == SAM_HASH
    assert census['total']['calls'] == 7241438
    assert census['validation']['sampled_calls'] == 36208
    assert census['validation']['cigar_consumption_beyond_input'] == 0
    assert census['total']['calls'] == sum(row['calls'] for row in census['by_kind'].values())
    close(census['total']['seconds'], sum(row['seconds'] for row in census['by_kind'].values()))
    print(f'Verified {len(paths)} imported artifacts, all bundle hashes, final trial aggregates, and 1x/10x/census receipt consistency.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--evidence', type=Path, default=REPO / 'docs/validation/2026-10-13-macos-arm64/evidence')
    verify(parser.parse_args().evidence)
