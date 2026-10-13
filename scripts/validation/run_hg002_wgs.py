#!/usr/bin/env python3
"""Acquire whole-read HG002 ONT prefixes; time and fully validate native sorted BAM."""
import argparse
import collections
import math
import platform
import urllib.request
from benchmark_sam import fastq_inventory, timing_metrics, time_command
import concurrent.futures
import datetime
import gzip
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import threading
import time
import traceback

REPO = Path(__file__).resolve().parents[2]
INPUTS = json.loads(Path(__file__).with_name('hg002_r1041_sources.json').read_text())
EXPECTED_REFERENCE_SHA = INPUTS['reference_sha256']
DENOMINATOR = INPUTS['primary_nuclear_bases']
WORK = OUT = DATA = BONEFISH = REFERENCE = OLD_BAM = None
REFERENCE_ROWS = []
SAMTOOLS = 'samtools'
COVERAGE, LABEL, TARGET, THREADS, SORT_MEM = 10, '10', 10 * DENOMINATOR, 8, '1G'
OLD_LENGTHS = {}
START = time.monotonic()
LOCK = threading.Lock()
ACTIVE = set()
STATUS = dict(stage='initializing', pid=os.getpid(), parts={})
COMP = bytes.maketrans(b'ACGTNacgtn', b'TGCANtgcan')
CIGAR_RE = re.compile(rb'(\d+)([MIDNSHP=X])')

def stamp():
    return datetime.datetime.now(datetime.timezone.utc).isoformat()

def save_json(path, value):
    temp = Path(str(path) + '.tmp')
    temp.write_text(json.dumps(value, indent=2) + '\n')
    temp.replace(path)

def emit(event, **fields):
    print(json.dumps(dict(at=stamp(), event=event, **fields)), flush=True)

def update(**fields):
    with LOCK:
        STATUS.update(fields)
        STATUS.update(updated_at_utc=stamp(), elapsed_seconds=round(time.monotonic() - START, 1))
        save_json(OUT / 'STATUS.json', STATUS)

def part_update(part, **fields):
    with LOCK:
        STATUS['parts'].setdefault(str(part), {}).update(fields)
        STATUS.update(updated_at_utc=stamp(), elapsed_seconds=round(time.monotonic() - START, 1))
        save_json(OUT / 'STATUS.json', STATUS)

def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda: f.read(8 << 20), b''):
            h.update(block)
    return h.hexdigest()

def reference_contigs(path):
    """Read names and lengths from FASTA, including ordinary gzip (not only BGZF)."""
    rows = []
    opener = gzip.open if path.suffix == '.gz' else open
    with opener(path, 'rb') as source:
        name, length = None, 0
        for line in source:
            if line.startswith(b'>'):
                if name is not None:
                    rows.append((name, length))
                name, length = line[1:].split()[0].decode(), 0
            else:
                if name is None and line.strip():
                    raise ValueError('Sequence before first FASTA header')
                length += len(line.strip())
        if name is not None:
            rows.append((name, length))
    if not rows or any(length <= 0 for _, length in rows) or len({name for name, _ in rows}) != len(rows):
        raise ValueError('Empty or duplicate FASTA contig')
    return rows

def read_fingerprint(seq, qual):
    h = hashlib.blake2b(digest_size=16)
    h.update(seq)
    h.update(b'\0')
    h.update(qual)
    return h.hexdigest()

def checked(args, **kwargs):
    r = subprocess.run([str(a) for a in args], capture_output=True, **kwargs)
    if r.returncode:
        error = r.stderr if isinstance(r.stderr, str) else r.stderr.decode(errors='replace')
        raise RuntimeError(f'Command returned {r.returncode}: {args}\n{error[-5000:]}')
    return r

def stop(proc):
    if proc.poll() is None:
        proc.terminate()
        try:
            proc.wait(timeout=10)
        except subprocess.TimeoutExpired:
            proc.kill()
            proc.wait()
    with LOCK:
        ACTIVE.discard(proc)

def acquire(part, source, target):
    fq = DATA / f'HG002_R1041_UL.part{part}.nominal{LABEL}x.fastq.gz'
    partial = Path(str(fq) + '.partial')
    receipt = OUT / f'part{part}.receipt.json'
    if receipt.exists():
        r = json.loads(receipt.read_text())
        assert r['source'] == source and r['target_bases'] == target
        assert digest(fq) == r['fastq_sha256']
        part_update(part, stage='complete', reads=r['reads'], bases=r['bases'])
        return r
    if partial.exists():
        # Preserve the last failed generated prefix before a deliberate retry.
        partial.rename(Path(str(partial) + f'.previous-{time.time_ns()}'))
    began = time.monotonic()
    count = bases = skipped = 0
    headers = []
    source_proc = None
    log_path = OUT / f'part{part}.remote.log'
    part_update(part, stage='retrieving', reads=0, bases=0, target_bases=target)
    try:
        with log_path.open('wb') as err, partial.open('wb') as raw:
            source_proc = subprocess.Popen([SAMTOOLS, 'view', '--no-PG', '-h', source['url']],
                                           stdout=subprocess.PIPE, stderr=err)
            with LOCK:
                ACTIVE.add(source_proc)
            with gzip.GzipFile(filename='', fileobj=raw, mode='wb', compresslevel=1, mtime=0) as out:
                for line in source_proc.stdout:
                    if line.startswith(b'@'):
                        headers.append(line.decode().rstrip())
                        continue
                    fields = line.rstrip(b'\n').split(b'\t', 11)
                    if len(fields) < 11:
                        raise RuntimeError(f'Malformed source SAM in part {part}')
                    flag = int(fields[1])
                    if flag & 0x900:
                        skipped += 1
                        continue
                    if not flag & 4 or flag & 16:
                        raise RuntimeError(f'Expected unaligned forward source record, part {part}, flag {flag}')
                    seq, qual = fields[9:11]
                    if seq == b'*' or qual == b'*' or not seq or len(seq) != len(qual):
                        raise RuntimeError(f'Invalid sequence/quality in part {part}')
                    out.write(b'@' + fields[0] + b'\n' + seq + b'\n+\n' + qual + b'\n')
                    count += 1
                    bases += len(seq)
                    if count % 2000 == 0:
                        part_update(part, reads=count, bases=bases, compressed_bytes=raw.tell(),
                                    seconds=round(time.monotonic() - began, 1))
                    if bases >= target:
                        break
            stop(source_proc)  # Intentional prefix termination; not a full-source EOF.
            source_proc.stdout.close()
        if bases < target:
            raise RuntimeError(f'Part {part} ended at {bases} bases, target {target}; {log_path.read_text()[-3000:]}')
        if not any(h.startswith('@RG\t') and '\tPL:ONT' in h for h in headers):
            raise RuntimeError(f'Part {part} is missing ONT read-group provenance')
        if not any('dorado' in h and '0.4.0' in h and 'sup@v4.1.0' in h for h in headers):
            raise RuntimeError(f'Part {part} does not match the expected Dorado model')
        partial.replace(fq)
        r = dict(part=part, source=source, source_headers=headers, target_bases=target,
                 reads=count, bases=bases, skipped_nonprimary=skipped,
                 fastq_path=str(fq), fastq_bytes=fq.stat().st_size, fastq_sha256=digest(fq),
                 seconds=round(time.monotonic() - began, 3), source_complete_file_downloaded=False,
                 sampling='First complete primary records in source file order. No length or quality filter.')
        save_json(receipt, r)
        part_update(part, stage='complete', reads=count, bases=bases, seconds=r['seconds'], compressed_bytes=r['fastq_bytes'])
        emit('part_complete', part=part, reads=count, bases=bases, seconds=r['seconds'])
        return r
    except BaseException as e:
        part_update(part, stage='failed', reads=count, bases=bases, error=str(e))
        raise
    finally:
        if source_proc is not None:
            stop(source_proc)

def prepare_dataset():
    update(stage='retrieving_reads')
    each, remainder = divmod(TARGET, 4)
    with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:
        futures = [pool.submit(acquire, i + 1, s, each + (i < remainder)) for i, s in enumerate(INPUTS['sources'])]
        parts = [f.result() for f in futures]
    update(stage='validating_input_fastqs')
    expected = {}
    lengths = []
    bases = 0
    for p in parts:
        part_reads = part_bases = 0
        with gzip.open(p['fastq_path'], 'rb') as f:
            while True:
                header = f.readline()
                if not header:
                    break
                seq = f.readline().rstrip(b'\r\n')
                plus = f.readline()
                qual = f.readline().rstrip(b'\r\n')
                if not header.startswith(b'@') or not plus.startswith(b'+') or not seq or len(seq) != len(qual):
                    raise RuntimeError(f'Invalid FASTQ framing in part {p["part"]}')
                if seq.translate(None, b'ACGTN') or any(q < 33 or q > 126 for q in qual):
                    raise RuntimeError('Unexpected source alphabet; sequence-preservation hash needs review')
                name = header[1:].split()[0].decode()
                if name in expected:
                    raise RuntimeError(f'Duplicate source read ID: {name}')
                expected[name] = [len(seq), read_fingerprint(seq, qual)]
                lengths.append(len(seq))
                part_reads += 1
                part_bases += len(seq)
                bases += len(seq)
                if part_reads % 10000 == 0:
                    update(validation_part=p['part'], validated_reads=len(expected), validated_bases=bases)
        if (part_reads, part_bases) != (p['reads'], p['bases']):
            raise RuntimeError('FASTQ counts do not match source receipt')
    assert bases >= TARGET
    assert expected, 'No input reads'
    old_names = OLD_LENGTHS
    assert set(old_names) <= set(expected), 'Expanded prefixes lost previously sampled reads'
    assert all(expected[k][0] == v for k, v in old_names.items())
    cumulative = 0
    n50 = None
    for length in sorted(lengths, reverse=True):
        cumulative += length
        if cumulative * 2 >= bases:
            n50 = length
            break
    provenance = dict(sample_label='HG002', platform='ONT', flowcell='R10.4.1',
                      basecaller='Dorado 0.4.0+0aaf16d', model='dna_r10.4.1_e8.2_400bps_sup@v4.1.0',
                      sampling='Equal read-base prefixes from four source BAMs, retaining whole primary records; file-order sampling, not uniform random sampling.',
                      target_bases=TARGET, reads=len(expected), bases=bases,
                      nominal_read_base_coverage=bases / DENOMINATOR,
                      coverage_denominator_bases=DENOMINATOR,
                      coverage_denominator_contigs=INPUTS['primary_nuclear_contigs'],
                      coverage_definition='Sequenced bases / full GRCh38 chr1-22+X+Y lengths including Ns; mapped depth is measured separately.',
                      read_length_min=min(lengths), read_length_max=max(lengths),
                      read_length_mean=bases / len(expected), read_length_n50=n50,
                      fastq_paths=[p['fastq_path'] for p in parts], source_parts=parts,
                      reference_path=str(REFERENCE), reference_sha256=digest(REFERENCE),
                      fastq_crc_and_record_validation_passed=True, duplicate_read_ids=0,
                      previous_subset_reads_retained=len(old_names), created_at_utc=stamp())
    assert provenance['reference_sha256'] == EXPECTED_REFERENCE_SHA
    save_json(OUT / 'DATASET_PROVENANCE.json', provenance)
    save_json(WORK / 'expected-read-fingerprints.json', expected)
    emit('dataset_ready', reads=len(expected), bases=bases, nominal_coverage=bases / DENOMINATOR, n50=n50)
    return provenance, expected, set(old_names)

def old_record_hashes(old_names):
    update(stage='reading_previous_subset_output_for_comparison')
    hashes = collections.defaultdict(list)
    with (OUT / 'previous-bam-decode.log').open('wb') as err:
        proc = subprocess.Popen([SAMTOOLS, 'view', OLD_BAM], stdout=subprocess.PIPE, stderr=err)
        try:
            for line in proc.stdout:
                name = line.split(b'\t', 1)[0].decode()
                assert name in old_names
                hashes[name].append(hashlib.sha256(line).hexdigest())
            if proc.wait():
                raise RuntimeError('Could not fully decode previous 1x BAM')
        finally:
            stop(proc)
    assert set(hashes) == old_names
    return {k: sorted(v) for k, v in hashes.items()}

def align(dataset):
    update(stage='aligning_and_sorting', alignment_started_utc=stamp())
    partial = DATA / f'bonefish.HG002.{LABEL}x.partial.bam'
    if partial.exists():
        raise RuntimeError('An existing partial BAM was retained; review before rerunning alignment')
    command = [str(BONEFISH), '-a', '-t', str(THREADS), '--sort', '--sort-mem', SORT_MEM,
               '--sort-tmp', str(DATA / f'bonefish.{LABEL}x.sort'), '-o', str(partial),
               str(REFERENCE), *dataset['fastq_paths']]
    env = os.environ.copy()
    for key in ('BONEFISH_CENSUS', 'BONEFISH_CENSUS_SEQ_EVERY', 'BONEFISH_CENSUS_SEQ_MIN'):
        env.pop(key, None)
    env['BONEFISH_SORT_TIMING'] = '1'
    log_path = OUT / 'bonefish.run.log'
    began = time.monotonic()
    emit('alignment_started', command=command)
    with log_path.open('wb') as log:
        proc = subprocess.Popen([*time_command(), *command], stdout=subprocess.DEVNULL, stderr=log, env=env)
        with LOCK:
            ACTIVE.add(proc)
        try:
            while proc.poll() is None:
                update(alignment_seconds=round(time.monotonic() - began, 1),
                       bam_bytes_written=partial.stat().st_size if partial.exists() else 0,
                       sort_run_files=len(list(DATA.glob(f'bonefish.{LABEL}x.sort.*.run'))))
                try:
                    proc.wait(timeout=10)
                except subprocess.TimeoutExpired:
                    pass
            code = proc.returncode
        finally:
            stop(proc)
    elapsed = time.monotonic() - began
    text = log_path.read_text()
    run = dict(command=command, environment={'BONEFISH_SORT_TIMING': '1'}, exit_code=code,
               wall_seconds_including_index_sort_bam_and_bai=elapsed, threads=THREADS, sort_mem=SORT_MEM,
               **timing_metrics(text),
               warnings=[l for l in text.splitlines() if re.search(r'warn|panic|error|failed', l, re.I)],
               binary_path=str(BONEFISH), binary_sha256=digest(BONEFISH), partial_bam_path=str(partial),
               started_utc=STATUS['alignment_started_utc'], completed_utc=stamp())
    save_json(OUT / 'bonefish.run.json', run)
    if code or run['warnings']:
        raise RuntimeError(f'Bonefish returned {code} or warned: {run["warnings"]}; see {log_path}')
    emit('alignment_complete', seconds=elapsed, peak_rss_bytes=run['peak_rss_bytes'])
    return partial, run

def validate_bam(partial, expected, old_hashes):
    update(stage='validating_full_bam', validation_records=0)
    checked([SAMTOOLS, 'quickcheck', '-v', partial])
    reference_rows = REFERENCE_ROWS or reference_contigs(REFERENCE)
    reference_lengths = {k: int(v) for k, v in reference_rows}
    tids = {k: i for i, (k, _) in enumerate(reference_rows)}
    seen = set()
    old_current = collections.defaultdict(list)
    counts = collections.Counter()
    contigs = collections.Counter()
    query_aligned = nuclear_match_bases = 0
    previous_key = (-1, -1, -1)
    headers = []
    with (OUT / 'bam-decode.log').open('wb') as err:
        proc = subprocess.Popen([SAMTOOLS, 'view', '-h', partial], stdout=subprocess.PIPE, stderr=err)
        with LOCK:
            ACTIVE.add(proc)
        try:
            for line in proc.stdout:
                if line.startswith(b'@'):
                    headers.append(line.decode().rstrip())
                    continue
                fields = line.rstrip(b'\n').split(b'\t')
                assert len(fields) >= 11
                name = fields[0].decode()
                flag = int(fields[1])
                assert name in expected, f'Unexpected output read: {name}'
                mapped = not bool(flag & 4)
                primary = not bool(flag & 0x900)
                ref = fields[2].decode()
                pos = int(fields[3])
                key = (tids[ref], pos, bool(flag & 16)) if mapped else (len(tids), 0, bool(flag & 16))
                assert key >= previous_key, f'Coordinate sort violation: {name}'
                previous_key = key
                counts['records'] += 1
                counts['secondary'] += bool(flag & 0x100)
                counts['supplementary'] += bool(flag & 0x800)
                seq, qual = fields[9:11]
                assert seq == b'*' or qual == b'*' or len(seq) == len(qual)
                ops = CIGAR_RE.findall(fields[5]) if mapped else []
                if mapped:
                    assert ops and b''.join(n + op for n, op in ops) == fields[5]
                    ops = [(int(n), op) for n, op in ops]
                    rspan = sum(n for n, op in ops if op in (b'M', b'D', b'N', b'=', b'X'))
                    qspan = sum(n for n, op in ops if op in (b'M', b'I', b'S', b'=', b'X'))
                    assert pos >= 1 and rspan > 0 and pos + rspan - 1 <= reference_lengths[ref]
                    assert seq == b'*' or qspan == len(seq), f'SEQ/CIGAR mismatch: {name}'
                if primary:
                    assert name not in seen, f'Duplicate primary: {name}'
                    seen.add(name)
                    counts['primary'] += 1
                    counts['mapped_primary'] += mapped
                    counts['unmapped_primary'] += not mapped
                    assert seq != b'*' and qual != b'*' and len(seq) == expected[name][0]
                    if flag & 16:
                        seq, qual = seq.translate(COMP)[::-1], qual[::-1]
                    assert read_fingerprint(seq, qual) == expected[name][1], f'Input sequence/quality changed: {name}'
                    if mapped:
                        contigs[ref] += 1
                        query_aligned += sum(n for n, op in ops if op in (b'M', b'I', b'=', b'X'))
                        if ref in INPUTS['primary_nuclear_contigs']:
                            nuclear_match_bases += sum(n for n, op in ops if op in (b'M', b'=', b'X'))
                if name in old_hashes:
                    old_current[name].append(hashlib.sha256(line).hexdigest())
                if counts['records'] % 10000 == 0:
                    update(validation_records=counts['records'], validation_primary=len(seen))
            if proc.wait():
                raise RuntimeError('Full BAM decode failed')
        finally:
            stop(proc)
    assert seen == set(expected), f'Missing primary reads: {len(set(expected) - seen)}'
    sq = [(next(f[3:] for f in h.split('\t') if f.startswith('SN:')),
           int(next(f[3:] for f in h.split('\t') if f.startswith('LN:')))) for h in headers if h.startswith('@SQ\t')]
    assert sq == [(k, int(v)) for k, v in reference_rows]
    assert any(h.startswith('@HD\t') and '\tSO:coordinate' in h for h in headers)
    differences = [k for k, v in old_hashes.items() if sorted(old_current.get(k, [])) != v]
    assert not differences, f'{len(differences)} previous 1x reads have different alignment records'
    result = dict(quickcheck_passed=True, full_decode_passed=True, coordinate_sort_passed=True,
                  exact_primary_read_inventory=True, missing_primary_reads=0, unexpected_read_names=0,
                  input_sequences_and_qualities_preserved=True, cigar_and_reference_bounds_passed=True,
                  reference_header_matches=True, counts=dict(counts), primary_contig_read_counts=dict(contigs),
                  mapped_primary_fraction=counts['mapped_primary'] / counts['primary'],
                  primary_aligned_query_bases=query_aligned,
                  primary_nuclear_match_mismatch_bases=nuclear_match_bases,
                  primary_nuclear_mean_depth_from_cigar=nuclear_match_bases / DENOMINATOR,
                  previous_subset_reads_compared=len(old_hashes),
                  previous_subset_alignment_records_compared=sum(map(len, old_hashes.values())),
                  previous_subset_record_sets_identical=True if old_hashes else None,
                  comparison_method='Decoded SAM record SHA256 multisets per read; headers excluded. Same-platform scaling check, not x86/ARM qualification.')
    save_json(OUT / 'BAM_VALIDATION.json', result)
    final = DATA / f'bonefish.HG002.{LABEL}x.sorted.bam'
    assert not final.exists() and not Path(str(final) + '.bai').exists()
    assert Path(str(partial) + '.bai').is_file()
    partial.rename(final)
    Path(str(partial) + '.bai').rename(Path(str(final) + '.bai'))
    emit('full_bam_validation_passed', counts=dict(counts), previous_subset_records_identical=True if old_hashes else None)
    return final, result

def coverage_and_receipt(bam):
    update(stage='measuring_coverage_and_checksums')
    coverage = checked([SAMTOOLS, 'coverage', '--ff', '0xF04', '-q', '0', '-Q', '0', bam], text=True).stdout
    (OUT / 'PRIMARY_COVERAGE.tsv').write_text(coverage)
    rows = []
    for line in coverage.splitlines():
        if not line or line.startswith('#'):
            continue
        f = line.split('\t')
        rows.append(dict(contig=f[0], length=int(f[2]) - int(f[1]) + 1, reads=int(f[3]),
                         covered_bases=int(f[4]), covered_percent=float(f[5]), mean_depth=float(f[6]),
                         mean_base_quality=float(f[7]), mean_mapping_quality=float(f[8])))
    nuclear = [r for r in rows if r['contig'] in INPUTS['primary_nuclear_contigs']]
    assert len(nuclear) == 24
    region = 'chr20:1-64444167'
    n_region = int(checked([SAMTOOLS, 'view', '-c', bam, region], text=True).stdout)
    assert n_region > 0
    return dict(bam_path=str(bam), bam_bytes=bam.stat().st_size, bam_sha256=digest(bam),
                bai_path=str(bam) + '.bai', bai_sha256=digest(Path(str(bam) + '.bai')),
                index_query=region, index_query_records=n_region,
                primary_only_exclude_flags='0xF04', minimum_mapping_quality=0, minimum_base_quality=0,
                nuclear_mean_depth=sum(r['mean_depth'] * r['length'] for r in nuclear) / DENOMINATOR,
                nuclear_covered_fraction=sum(r['covered_bases'] for r in nuclear) / DENOMINATOR,
                per_contig_primary_coverage=rows)

def source_headers():
    sources = []
    for source in INPUTS['sources']:
        request = urllib.request.Request(source['url'], method='HEAD')
        with urllib.request.urlopen(request, timeout=30) as response:
            actual = dict(url=source['url'], status=response.status,
                          bytes=int(response.headers['Content-Length']),
                          etag=response.headers.get('ETag'),
                          last_modified=response.headers.get('Last-Modified'),
                          accept_ranges=response.headers.get('Accept-Ranges'))
        actual['matches_previous'] = all(actual[k] == source[k] for k in ('bytes', 'etag', 'last_modified'))
        if not actual['matches_previous']:
            raise RuntimeError(f'Public source changed since the recorded benchmark: {source["url"]}')
        sources.append(actual)
    return sources


def main():
    global WORK, OUT, DATA, BONEFISH, REFERENCE, OLD_BAM, OLD_LENGTHS, SAMTOOLS
    global COVERAGE, LABEL, TARGET, THREADS, SORT_MEM
    global REFERENCE_ROWS
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--bonefish', type=Path, required=True)
    ap.add_argument('--reference', type=Path, required=True)
    ap.add_argument('--data-dir', type=Path, required=True)
    ap.add_argument('--output-dir', type=Path, required=True)
    ap.add_argument('--work-dir', type=Path, required=True)
    ap.add_argument('--samtools', default='samtools')
    ap.add_argument('--coverage', type=float, default=10)
    ap.add_argument('--threads', type=int, default=8)
    ap.add_argument('--sort-mem', default='1G')
    ap.add_argument('--minimum-free-gib', type=float, help='Default: 15 GiB per nominal x; an estimate, not a capacity guarantee')
    ap.add_argument('--expected-binary-sha256', help='Optional exact build check before running')
    ap.add_argument('--previous-bam', type=Path, help='Optional earlier sorted output for exact per-read record comparison')
    ap.add_argument('--previous-reads', type=Path, nargs='+', help='FASTQs used for --previous-bam (supply both)')
    args = ap.parse_args()
    if not math.isfinite(args.coverage) or args.coverage <= 0 or args.threads < 1:
        ap.error('coverage must be finite and positive; threads must be positive')
    if bool(args.previous_bam) != bool(args.previous_reads):
        ap.error('Supply --previous-bam and --previous-reads together')
    min_free = args.minimum_free_gib if args.minimum_free_gib is not None else 15 * args.coverage
    if not math.isfinite(min_free) or min_free < 0:
        ap.error('minimum-free-gib must be finite and nonnegative')
    WORK, OUT, DATA = (p.resolve() for p in (args.work_dir, args.output_dir, args.data_dir))
    BONEFISH, REFERENCE = args.bonefish.resolve(), args.reference.resolve()
    OLD_BAM = args.previous_bam.resolve() if args.previous_bam else None
    SAMTOOLS, COVERAGE, THREADS, SORT_MEM = args.samtools, args.coverage, args.threads, args.sort_mem
    LABEL, TARGET = format(COVERAGE, 'g'), max(1, round(COVERAGE * DENOMINATOR))
    for directory in (WORK, OUT, DATA):
        directory.mkdir(parents=True, exist_ok=True)
    if (OUT / 'RESULTS.json').exists() or any(DATA.glob(f'bonefish.HG002.{LABEL}x.*.bam*')):
        raise RuntimeError('Existing alignment/results retained; use a new data/output directory')
    if shutil.disk_usage(DATA).free < min_free * (1 << 30):
        raise RuntimeError(f'Less than {min_free:g} GiB free in data directory')
    binary_sha = digest(BONEFISH)
    if args.expected_binary_sha256 and binary_sha != args.expected_binary_sha256:
        raise RuntimeError('Binary SHA256 does not match the requested build')
    if digest(REFERENCE) != EXPECTED_REFERENCE_SHA:
        raise RuntimeError('Expected the recorded GRCh38 no-alt reference (compressed-file SHA256 mismatch)')
    # Read bounds outside timing; samtools faidx cannot index ordinary gzip references.
    REFERENCE_ROWS = reference_contigs(REFERENCE)
    assert sum(n for _, n in REFERENCE_ROWS) == INPUTS['all_reference_bases']
    assert sum(n for name, n in REFERENCE_ROWS if name in INPUTS['primary_nuclear_contigs']) == DENOMINATOR
    if args.previous_reads:
        OLD_LENGTHS = {name.decode(): value[0] for name, value in fastq_inventory(args.previous_reads).items()}
    save_json(OUT / 'SOURCE_HTTP_HEADERS.json', source_headers())
    update(stage='preflight_complete', target_bases=TARGET, target_nominal_coverage=COVERAGE,
           data_directory=str(DATA), free_data_bytes=shutil.disk_usage(DATA).free)
    emit('job_started', target_bases=TARGET, data_directory=str(DATA))
    dataset, expected, old_names = prepare_dataset()
    old_hashes = old_record_hashes(old_names) if OLD_BAM else {}
    partial, run = align(dataset)
    bam, validation = validate_bam(partial, expected, old_hashes)
    bam_receipt = coverage_and_receipt(bam)
    duration = run['wall_seconds_including_index_sort_bam_and_bai']
    git = checked(['git', '-C', REPO, 'rev-parse', 'HEAD'], text=True).stdout.strip() if (REPO / '.git').exists() else None
    result = dict(dataset=dataset, bonefish_run=run, validation=validation, bam_receipt=bam_receipt,
                  throughput_reads_per_second=dataset['reads'] / duration,
                  throughput_sequenced_bases_per_second=dataset['bases'] / duration,
                  machine=platform.platform(), samtools_version=checked([SAMTOOLS, '--version-only'], text=True).stdout.strip(),
                  runner_source_head=git, all_checks_passed=True, completed_utc=stamp(),
                  total_job_seconds=time.monotonic() - START,
                  scope='One HG002 ONT file-order prefix sample. Fresh index, mapping, native coordinate sort, BAM compression and BAI are timed. Acquisition, input/output checks, coverage and hashes are outside timing. Structural and preservation checks do not establish accuracy or cross-platform equality.')
    save_json(OUT / 'RESULTS.json', result)
    checksums = ''.join(f'{digest(p)}  {p.name}\n' for p in sorted(OUT.iterdir())
                        if p.is_file() and p.name not in ('SHA256SUMS.txt', 'STATUS.json'))
    (OUT / 'SHA256SUMS.txt').write_text(checksums)
    update(stage='complete', result_path=str(OUT / 'RESULTS.json'))
    emit('job_complete', reads=dataset['reads'], bases=dataset['bases'], seconds=duration,
         mean_depth=bam_receipt['nuclear_mean_depth'], covered_fraction=bam_receipt['nuclear_covered_fraction'])


if __name__ == '__main__':
    try:
        main()
    except SystemExit:
        raise
    except BaseException as error:
        if OUT is not None and OUT.exists():
            update(stage='failed', error=str(error), traceback=traceback.format_exc())
            save_json(OUT / 'FAILURE.json', dict(at=stamp(), error=str(error), traceback=traceback.format_exc()))
        with LOCK:
            processes = list(ACTIVE)
        for proc in processes:
            stop(proc)
        emit('job_failed', error=str(error))
        raise
