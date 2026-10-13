#!/usr/bin/env python3
"""Time one SAM run; check every primary read and hash output after timing."""
import argparse
import datetime
import gzip
import hashlib
import json
from pathlib import Path
import platform
import re
import subprocess
import time


def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda: f.read(8 << 20), b''):
            h.update(block)
    return h.hexdigest()


def read_fingerprint(seq, qual):
    h = hashlib.blake2b(digest_size=16)
    h.update(seq)
    h.update(b'\0')
    h.update(qual)
    return h.hexdigest()


def fastq_inventory(paths):
    expected = {}
    for path in paths:
        opener = gzip.open if path.suffix == '.gz' else open
        with opener(path, 'rb') as f:
            while header := f.readline():
                seq, plus, qual = f.readline().rstrip(b'\r\n'), f.readline(), f.readline().rstrip(b'\r\n')
                if not header.startswith(b'@') or not plus.startswith(b'+') or not seq or len(seq) != len(qual):
                    raise ValueError(f'Invalid FASTQ framing in {path}')
                if seq.translate(None, b'ACGTN') or any(q < 33 or q > 126 for q in qual):
                    raise ValueError(f'Invalid FASTQ sequence/quality in {path}')
                name = header[1:].split()[0]
                if name in expected:
                    raise ValueError(f'Duplicate FASTQ read ID: {name!r}')
                expected[name] = (len(seq), read_fingerprint(seq, qual))
    if not expected:
        raise ValueError('No input reads')
    return expected


def time_command():
    # BSD time reports bytes; GNU time reports maximum RSS in KiB.
    return ['/usr/bin/time', '-l' if platform.system() == 'Darwin' else '-v']


def timing_metrics(log):
    rss = re.search(r'^\s*(\d+)\s+maximum resident set size', log, re.M)
    if rss:
        rss_bytes = int(rss[1])
    else:
        rss = re.search(r'Maximum resident set size \(kbytes\):\s*(\d+)', log)
        rss_bytes = int(rss[1]) * 1024 if rss else None
    result = {'peak_rss_bytes': rss_bytes}
    for key, label in [('instructions_retired', 'instructions retired'),
                       ('peak_memory_footprint_bytes', 'peak memory footprint')]:
        match = re.search(r'^\s*(\d+)\s+' + label, log, re.M)
        result[key] = int(match[1]) if match else None
    stages = {}
    for label in ['Sequence Sketching', 'Seeding (Lookup)', 'Chaining', 'Alignment (Ext)', 'Post-chain', 'Total Measured']:
        match = re.search(re.escape(label) + r':\s+([\d.]+)s', log)
        if match:
            stages[label] = float(match[1])
    result['worker_stage_sums_seconds'] = stages
    return result


def validate_sam(path, expected):
    h = hashlib.sha256()
    seen = set()
    records = mapped = 0
    comp = bytes.maketrans(b'ACGTN', b'TGCAN')
    with path.open('rb') as f:
        for line in f:
            if line.startswith(b'@PG'):
                continue
            h.update(line)
            if line.startswith(b'@'):
                continue
            fields = line.rstrip(b'\n').split(b'\t')
            if len(fields) < 11 or fields[0] not in expected:
                raise ValueError('Malformed SAM record or unexpected read ID')
            records += 1
            flag = int(fields[1])
            if not flag & 0x900:
                name, seq, qual = fields[0], fields[9], fields[10]
                if name in seen or seq == b'*' or qual == b'*' or len(seq) != len(qual):
                    raise ValueError('Duplicate primary or invalid primary sequence/quality')
                if flag & 16:
                    seq, qual = seq.translate(comp)[::-1], qual[::-1]
                if (len(seq), read_fingerprint(seq, qual)) != expected[name]:
                    raise ValueError(f'Primary sequence/quality differs from input: {name!r}')
                seen.add(name)
                mapped += not bool(flag & 4)
    if seen != set(expected):
        raise ValueError(f'Missing primary reads: {len(set(expected) - seen)}')
    return dict(sam_sha256_without_pg=h.hexdigest(), records=records, primary=len(seen),
                unique_primary_names=len(seen), mapped_primary=mapped,
                missing_primary_names=0, unexpected_primary_names=0,
                input_sequences_and_qualities_preserved=True)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('binary', type=Path)
    p.add_argument('--reference', type=Path, required=True)
    p.add_argument('--reads', type=Path, nargs='+', required=True)
    p.add_argument('--output-dir', type=Path, required=True)
    p.add_argument('--variant', choices=['control', 'candidate'], required=True)
    p.add_argument('--trial', type=int, required=True)
    p.add_argument('--threads', type=int, default=8)
    p.add_argument('--expected-sam-sha256', help='Use the original control hash to enforce output preservation')
    p.add_argument('--census-prefix', type=Path, help='For a census-feature build; sets BONEFISH_CENSUS and SEQ_EVERY=200')
    args = p.parse_args()
    if args.threads < 1 or args.trial < 1:
        p.error('threads and trial must be positive')
    args.output_dir.mkdir(parents=True, exist_ok=True)
    stem = args.output_dir / f'{args.variant}-{args.trial}'
    sam, log_path, result_path = (Path(str(stem) + ext) for ext in ('.sam', '.log', '.json'))
    if any(path.exists() for path in (sam, log_path, result_path)):
        p.error('Trial output exists; use a new trial number/directory')
    expected = fastq_inventory(args.reads)
    binary_hash = digest(args.binary)
    inputs = {str(path): digest(path) for path in [args.reference, *args.reads]}
    command = [str(args.binary.resolve()), '-a', '-t', str(args.threads),
               str(args.reference.resolve()), *(str(path.resolve()) for path in args.reads)]
    import os
    env = os.environ.copy()
    # Keep instrumentation explicitly controlled, even if the invoking shell sets it.
    env.pop('BONEFISH_CENSUS', None)
    env.pop('BONEFISH_CENSUS_SEQ_EVERY', None)
    env.pop('BONEFISH_CENSUS_SEQ_MIN', None)
    if args.census_prefix:
        if any(Path(str(args.census_prefix) + ext).exists() for ext in ('.calls.bin', '.seqs.bin')):
            p.error('Census files exist; use a new prefix')
        args.census_prefix.parent.mkdir(parents=True, exist_ok=True)
        env.update(BONEFISH_CENSUS=str(args.census_prefix.resolve()), BONEFISH_CENSUS_SEQ_EVERY='200')
    started = datetime.datetime.now(datetime.timezone.utc).isoformat()
    with sam.open('xb') as out, log_path.open('xb') as log:
        t0 = time.perf_counter()
        run = subprocess.run([*time_command(), *command], stdout=out, stderr=log, env=env)
        elapsed = time.perf_counter() - t0
    text = log_path.read_text()
    result = dict(run=f'{args.variant}-{args.trial}', started_utc=started, command=command,
                  binary_sha256=binary_hash, input_sha256=inputs, returncode=run.returncode,
                  threads=args.threads, elapsed_seconds_including_index_and_sam=elapsed,
                  census_instrumented=bool(args.census_prefix),
                  timing_scope='Fresh FASTA-built index, mapping and SAM output. Input/output hashing and read validation are outside the timer.',
                  warnings=[line for line in text.splitlines() if re.search(r'warn|panic|error|failed', line, re.I)],
                  **timing_metrics(text))
    try:
        if run.returncode:
            raise ValueError(f'Alignment returned {run.returncode}; see {log_path}')
        result.update(validate_sam(sam, expected))
        result['same_as_original_release_sam'] = (result['sam_sha256_without_pg'] == args.expected_sam_sha256
                                                  if args.expected_sam_sha256 else None)
        if args.expected_sam_sha256 and not result['same_as_original_release_sam']:
            raise ValueError('Normalized SAM differs from the expected original control hash')
        result['all_checks_passed'] = True
    except ValueError as error:
        result.update(all_checks_passed=False, error=str(error))
        raise
    finally:
        result_path.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
