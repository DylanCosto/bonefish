from pathlib import Path
import collections
import datetime
import gzip
import hashlib
import json
import re
import shutil

SESSION = Path('${SESSION}')
OUT = SESSION / 'outputs/bonefish-hg001-parity'
WORK = SESSION / 'work/bonefish-hg001-parity'
RUN = WORK / 'default-sam/candidate-1.json'
SAM = RUN.with_suffix('.sam')
PREFIX = 'deed7949d35f'
CIGAR = re.compile(rb'(\d+)([MIDNSHP=X])')
run = json.loads(RUN.read_text())
assert run['all_checks_passed'] and run['returncode'] == 0 and not run['warnings']
download = json.loads((OUT / 'DOWNLOAD_STATUS.json').read_text())
reference = json.loads((OUT / 'REFERENCE_RECEIPT.json').read_text())
assert run['input_sha256'][download['local_file']] == download['hashes']['sha256']
assert run['input_sha256'][reference['internal_copy']] == reference['sha256']

# Independently read the FASTA to qualify all SAM @SQ names, order and lengths.
ref_rows = []
with gzip.open(reference['internal_copy'], 'rb') as source:
    name, length = None, 0
    for line in source:
        if line.startswith(b'>'):
            if name is not None:
                ref_rows.append((name, length))
            name, length = line[1:].split()[0], 0
        else:
            length += len(line.strip())
    if name is not None:
        ref_rows.append((name, length))
assert len(ref_rows) == 195 and sum(n for _, n in ref_rows) == 3099922541
references = dict(ref_rows)
sam_rows = []
counts = collections.Counter()
lengths = []
hashes = {key: hashlib.sha256() for key in ['raw_sam', 'sam_without_pg', 'records_without_all_headers']}

def check_group(records):
    primary = [fields for fields in records if not int(fields[1]) & 0x900]
    assert len(primary) == 1
    query_length = len(primary[0][9])
    lengths.append(query_length)
    for fields in records:
        counts['records'] += 1
        flag = int(fields[1])
        counts['secondary' if flag & 0x100 else 'supplementary' if flag & 0x800 else 'primary'] += 1
        if not flag & 0x900:
            counts['unmapped_primary' if flag & 4 else 'mapped_primary'] += 1
        seq, qual = fields[9:11]
        assert seq == b'*' or qual == b'*' or len(seq) == len(qual)
        if flag & 4:
            continue
        ops = CIGAR.findall(fields[5])
        assert ops and b''.join(n + op for n, op in ops) == fields[5]
        ops = [(int(n), op) for n, op in ops]
        assert all(n > 0 for n, _ in ops)
        qspan = sum(n for n, op in ops if op in (b'M', b'I', b'S', b'=', b'X'))
        hard = sum(n for n, op in ops if op == b'H')
        rspan = sum(n for n, op in ops if op in (b'M', b'D', b'N', b'=', b'X'))
        assert qspan + hard == query_length, 'Full-query CIGAR consumption mismatch'
        assert seq == b'*' or qspan == len(seq), 'CIGAR/SEQ length mismatch'
        position = int(fields[3])
        assert fields[2] in references and position >= 1 and rspan > 0
        assert position + rspan - 1 <= references[fields[2]], 'Reference bounds exceeded'

records, name = [], None
with SAM.open('rb') as source:
    for line in source:
        hashes['raw_sam'].update(line)
        if not line.startswith(b'@PG'):
            hashes['sam_without_pg'].update(line)
        if line.startswith(b'@'):
            if line.startswith(b'@SQ\t'):
                fields = dict(f.split(b':', 1) for f in line.rstrip(b'\n').split(b'\t')[1:])
                sam_rows.append((fields[b'SN'], int(fields[b'LN'])))
            continue
        hashes['records_without_all_headers'].update(line)
        fields = line.rstrip(b'\n').split(b'\t')
        assert len(fields) >= 11
        if name is not None and fields[0] != name:
            check_group(records)
            records = []
        name = fields[0]
        records.append(fields)
    if records:
        check_group(records)
assert sam_rows == ref_rows
assert counts['records'] == run['records'] and counts['primary'] == run['primary']
assert counts['mapped_primary'] == run['mapped_primary']
assert counts['mapped_primary'] + counts['unmapped_primary'] == run['primary']
digests = {name: h.hexdigest() for name, h in hashes.items()}
assert digests['sam_without_pg'] == run['sam_sha256_without_pg']
total_bases = sum(lengths)
cumulative = 0
for length in sorted(lengths, reverse=True):
    cumulative += length
    if cumulative >= (total_bases + 1) // 2:
        n50 = length
        break
matches = [name for name, digest in digests.items() if digest.startswith(PREFIX)]
report = {
    'status': 'complete',
    'completed_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(),
    'tested_source_commit': 'c1b7eb074b49806a86bbd0e48c617957ac436618',
    'branch_documentation_commit': '22619975719bbf59ea816e6b6349bb8ad5fee170',
    'run': run,
    'dataset': {'drive_file_id': download['drive_file_id'], 'compressed_bytes': download['bytes'],
                'compressed_sha256': download['hashes']['sha256'], 'reads': len(lengths),
                'bases': total_bases, 'read_n50': n50,
                'identity_source': 'User-selected HG001 genome-wide sample and Drive filename; no independent donor or chemistry qualification.'},
    'counts': dict(counts),
    'sam_bytes': SAM.stat().st_size,
    'sam_sha256': digests,
    'validation': {'gzip_crc_and_fastq_framing_passed': True, 'exact_primary_inventory_passed': True,
                   'input_sequence_quality_fingerprints_passed': True, 'positive_complete_cigar_operations_passed': True,
                   'full_query_cigar_consumption_passed': True, 'cigar_sequence_lengths_passed': True,
                   'reference_header_matches_fasta': True, 'reference_bounds_passed': True},
    'expected_x86_hash_prefix_from_user': PREFIX,
    'expected_full_x86_hash': None,
    'prefix_matches': matches,
    'full_x86_arm_parity_established': False,
    'qualification': 'Only a 12-hex-digit expected x86 prefix was supplied. Exact x86 invocation, full hash, normalization and reference receipt were not provided. Matching a prefix is reported separately from full-SHA256 confirmation; a mismatch cannot be attributed to the kernel without those identities.',
    'storage': 'FASTQ, reference and SAM on internal APFS SSD; download, hashing and all validation outside the aligner timer.',
}
(OUT / 'PARITY_RESULTS.json').write_text(json.dumps(report, indent=2) + '\n')
for source in [RUN, RUN.with_suffix('.log'), WORK / 'default-sam.driver.log', WORK / 'download.log']:
    shutil.copy2(source, OUT / source.name)
summary = f'''# HG001 genome-wide Mac parity run

Tested source: `c1b7eb074b49806a86bbd0e48c617957ac436618` (clean NEON DP).
Drive FASTQ: `{download['drive_file_id']}`; compressed SHA256 `{download['hashes']['sha256']}`.
The download matched Drive's SHA256, SHA1 and MD5; file size was {download['bytes']:,} bytes.

The standard `bonefish -a -t 8 REFERENCE READS` run completed in
**{run['elapsed_seconds_including_index_and_sam']:.3f} s**, including fresh index and SAM writing.
FASTQ, reference and SAM were on the internal APFS SSD. This is one run, not a paired speed comparison.
Exact commands and input/binary hashes are in `candidate-1.json`.

Input: **{len(lengths):,} reads / {total_bases:,} bases**, read N50 {n50:,} bp.
Output: **{counts['records']:,} records**, {counts['mapped_primary']:,} mapped and
{counts['unmapped_primary']:,} unmapped primaries. Every input primary, sequence and quality
was preserved after restoring orientation. Full CIGAR consumption, SEQ lengths,
reference headers and bounds passed for every alignment record.

| Hash normalization | Complete SHA256 |
|---|---|
| Raw SAM, including all headers | `{digests['raw_sam']}` |
| SAM excluding only @PG | `{digests['sam_without_pg']}` |
| Exact record stream excluding all headers | `{digests['records_without_all_headers']}` |

Expected x86 prefix supplied: **`{PREFIX}`**.
Matching normalizations: **{', '.join(matches) if matches else 'none'}**.
The full x86 hash, exact invocation, normalization and reference identity remain unprovided;
full-SHA256 ARM/x86 parity is **not established** by this receipt.
The record-stream normalization excludes headers only; it preserves all records,
order, CIGARs, sequences, qualities and tags. No fields were removed to force a match.

The complete SAM is retained in internal scratch. Large FASTQ/reference/SAM files stay outside Git.
This sample alone does not reproduce the separate GIAB variant-calling benchmarks.
'''
(OUT / 'README.md').write_text(summary)
checksum_lines = [f'{hashlib.sha256(p.read_bytes()).hexdigest()}  {p.name}'
                  for p in sorted(OUT.iterdir()) if p.is_file() and p.name != 'SHA256SUMS.txt']
(OUT / 'SHA256SUMS.txt').write_text('\n'.join(checksum_lines) + '\n')
print(json.dumps({'elapsed_seconds': run['elapsed_seconds_including_index_and_sam'],
                  'reads': len(lengths), 'bases': total_bases, 'counts': dict(counts),
                  'hashes': digests, 'prefix_matches': matches}, indent=2), flush=True)
