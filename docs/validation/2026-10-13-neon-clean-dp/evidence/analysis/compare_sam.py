import collections
import hashlib
import itertools
import json
from pathlib import Path
import re

ROOT = Path('${SESSION}')
OUT = ROOT / 'outputs/bonefish-neon-clean-dp'
SCRATCH = ROOT / 'work/bonefish-neon-clean-dp/benchmark'
CIGAR = re.compile(rb'(\d+)([MIDNSHP=X])')

def reference_header(path):
    rows = []
    with path.open('rb') as source:
        for line in source:
            if not line.startswith(b'@'):
                break
            if line.startswith(b'@SQ\t'):
                fields = dict(field.split(b':', 1) for field in line.rstrip(b'\n').split(b'\t')[1:])
                rows.append((fields[b'SN'], int(fields[b'LN'])))
    assert len(rows) == 195 and sum(length for _, length in rows) == 3099922541
    return rows

def check_records(records, query_length, references):
    for line in records:
        fields = line.rstrip(b'\n').split(b'\t')
        flag = int(fields[1])
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
        assert qspan + hard == query_length, 'CIGAR consumption differs from full input query'
        assert seq == b'*' or qspan == len(seq), 'SEQ/CIGAR length mismatch'
        position = int(fields[3])
        assert fields[2] in references and position >= 1 and rspan > 0
        assert position + rspan - 1 <= references[fields[2]], 'Reference bounds exceeded'

def groups(path):
    name = None
    records = []
    with path.open('rb') as source:
        for line in source:
            if line.startswith(b'@'):
                continue
            current = line.split(b'\t', 1)[0]
            if name is not None and current != name:
                yield name, records
                records = []
            name = current
            records.append(line)
    if name is not None:
        yield name, records

def primary(records):
    found = []
    for line in records:
        fields = line.rstrip(b'\n').split(b'\t')
        if not int(fields[1]) & 0x900:
            found.append(fields)
    assert len(found) == 1
    return found[0]

def summary(fields):
    return dict(flag=int(fields[1]), contig=fields[2].decode(), position=int(fields[3]),
                mapq=int(fields[4]), cigar=fields[5].decode(),
                tags=[tag.decode() for tag in fields[11:] if tag.startswith((b'AS:', b'NM:', b'ms:', b'tp:'))])

for name in ['control-1', 'candidate-1']:
    assert json.loads((OUT / f'{name}.json').read_text())['all_checks_passed'], 'Compare only complete, validated SAMs'
ref_rows = reference_header(SCRATCH / 'control-1.sam')
assert reference_header(SCRATCH / 'candidate-1.sam') == ref_rows
references = dict(ref_rows)
counts = collections.Counter()
examples = []
material = []
with (OUT / 'CHANGED_READS.jsonl').open('w') as changes:
    pairs = itertools.zip_longest(groups(SCRATCH / 'control-1.sam'), groups(SCRATCH / 'candidate-1.sam'))
    for left, right in pairs:
        assert left is not None and right is not None
        name, old = left
        new_name, new = right
        assert name == new_name, 'Read-group order differs; use a keyed comparison instead'
        counts['reads_compared'] += 1
        counts['control_records'] += len(old)
        counts['candidate_records'] += len(new)
        old_set = sorted(hashlib.sha256(line).hexdigest() for line in old)
        new_set = sorted(hashlib.sha256(line).hexdigest() for line in new)
        equal = old_set == new_set
        counts['same_record_multiset_reads' if equal else 'changed_record_multiset_reads'] += 1
        if len(old) != len(new):
            counts['record_count_changed_reads'] += 1
        a, b = primary(old), primary(new)
        assert len(a[9]) == len(b[9])
        check_records(old, len(a[9]), references)
        check_records(new, len(b[9]), references)
        fa, fb = int(a[1]), int(b[1])
        ma, mb = not bool(fa & 4), not bool(fb & 4)
        for label, different in [('primary_record_changed', a != b), ('primary_flag_changed', a[1] != b[1]),
                                 ('primary_contig_changed', a[2] != b[2]), ('primary_position_changed', a[3] != b[3]),
                                 ('primary_strand_changed', bool(fa & 16) != bool(fb & 16)),
                                 ('primary_mapq_changed', a[4] != b[4]), ('primary_cigar_changed', a[5] != b[5]),
                                 ('primary_tags_changed', a[11:] != b[11:])]:
            counts[label] += different
        category = 'both_mapped' if ma and mb else 'both_unmapped' if not ma and not mb else 'gained_mapping' if mb else 'lost_mapping'
        counts[category] += 1
        different_placement = ma and mb and (a[2] != b[2] or bool(fa & 16) != bool(fb & 16))
        counts['both_mapped_primary_contig_or_strand_changed'] += different_placement
        if not equal:
            value = dict(read_id=name.decode(), control_records=len(old), candidate_records=len(new),
                         control_primary=summary(a), candidate_primary=summary(b))
            changes.write(json.dumps(value) + '\n')
            if len(examples) < 10:
                examples.append(value)
            if (different_placement or ma != mb) and len(material) < 20:
                material.append(value)
assert counts['reads_compared'] == 86553
report = dict(counts=dict(counts), example_changed_reads=examples, example_placement_changes=material,
              validation=dict(reference_headers_equal=True, complete_record_cigar_consumption_passed=True,
                              reference_bounds_passed=True, sequence_cigar_lengths_passed=True),
              method='Exact decoded SAM-line SHA256 multisets per read, excluding all headers; primary field comparisons. Both full SAMs independently passed exact input inventory and orientation-restored sequence/quality checks.',
              scope='Old ARM legacy/static DP versus this ARM clean/adaptive DP branch on one Mac. Changes do not establish accuracy or whole-WGS x86 parity; native kernel fixtures test the common clean semantics.')
(OUT / 'OUTPUT_COMPARISON.json').write_text(json.dumps(report, indent=2) + '\n')
print(json.dumps(report['counts'], indent=2), flush=True)
