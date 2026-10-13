#!/usr/bin/env python3
"""Decode Bonefish's 68-byte helper census and verify sampled input framing."""
import argparse
import collections
import hashlib
import heapq
import json
import mmap
import struct
from pathlib import Path

REC = struct.Struct('<BBH4i2I2i8I')
assert REC.size == 68
KINDS = ['left_extension', 'gap_fill', 'gap_rerun', 'right_extension']
BOUNDS = [32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384]
LABELS = ['0-31', '32-63', '64-127', '128-255', '256-511', '512-1023',
          '1024-2047', '2048-4095', '4096-8191', '8192-16383', '16384+']

def bucket(n):
    for i, bound in enumerate(BOUNDS):
        if n < bound:
            return i
    return len(BOUNDS)

def blank():
    return [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

def add(s, r):
    s[0] += 1
    s[1] += r[18]
    s[2] += r[7]
    s[3] += r[8]
    s[4] += bool(r[1] & 1)
    s[5] += bool(r[1] & 2)
    s[6] += bool(r[1] & 4)
    s[7] += bool(r[1] & 16)
    s[8] += r[18] == 0xffffffff
    s[9] = max(s[9], r[18])

def present(s, total_ns):
    return dict(calls=s[0], seconds=s[1] / 1e9,
                percent_of_recorded_helper_time=100 * s[1] / total_ns,
                mean_microseconds=s[1] / s[0] / 1000,
                mean_query_bases=s[2] / s[0], mean_target_bases=s[3] / s[0],
                zdropped=s[4], reach_end=s[5], sampled=s[6], guard_skipped=s[7],
                duration_saturated=s[8], max_milliseconds=s[9] / 1e6)

def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda: f.read(4 << 20), b''):
            h.update(block)
    return h.hexdigest()

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('prefix', type=Path)
    ap.add_argument('output', type=Path)
    a = ap.parse_args()
    calls = Path(str(a.prefix) + '.calls.bin')
    seqs = Path(str(a.prefix) + '.seqs.bin')
    total = blank()
    by_kind = [blank() for _ in KINDS]
    by_size = [[blank() for _ in LABELS] for _ in KINDS]
    by_flags = collections.defaultdict(blank)
    slow = []
    invalid_calls = 0
    with calls.open('rb') as f, mmap.mmap(f.fileno(), 0, access=mmap.ACCESS_READ) as mm:
        assert len(mm) % REC.size == 0, 'Incomplete call record'
        for index, r in enumerate(REC.iter_unpack(mm)):
            assert r[0] < 4 and r[2] == 0 and r[1] < 32, ('Invalid call header', index)
            add(total, r)
            add(by_kind[r[0]], r)
            add(by_size[r[0]][bucket(max(r[7], r[8]))], r)
            add(by_flags[(r[0], r[3])], r)
            if r[11] > r[7] or r[12] > r[8]:
                invalid_calls += 1
            candidate = (r[18], index, r)
            if len(slow) < 20:
                heapq.heappush(slow, candidate)
            elif candidate > slow[0]:
                heapq.heapreplace(slow, candidate)

    sample_count = 0
    sample_kinds = collections.Counter()
    with seqs.open('rb') as f, mmap.mmap(f.fileno(), 0, access=mmap.ACCESS_READ) as mm:
        off = 0
        while off < len(mm):
            assert off + 117 <= len(mm), ('Incomplete sample header', off)
            call_index, seq_id = struct.unpack_from('<IQ', mm, off)
            r = REC.unpack_from(mm, off + 12)
            qlen, tlen = struct.unpack_from('<II', mm, off + 109)
            assert r[0] < 4 and r[2] == 0 and r[1] & 4
            assert (qlen, tlen) == (r[7], r[8])
            end = off + 117 + qlen + tlen
            assert end <= len(mm), ('Truncated sample input', off)
            assert max(mm[off + 117:end], default=0) <= 4, 'Invalid nt4 input'
            sample_count += 1
            sample_kinds[KINDS[r[0]]] += 1
            off = end
        assert off == len(mm)
    assert sample_count == total[6], (sample_count, total[6])
    assert sample_count == (total[0] + 199) // 200
    assert sum(sample_kinds.values()) == sample_count
    assert all(sample_kinds[k] == s[6] for k, s in zip(KINDS, by_kind))
    assert invalid_calls == 0, 'Recorded CIGAR consumption exceeds input bounds'

    report = dict(
        format='Bonefish helper census, 68-byte little-endian records',
        input_files={p.name: dict(bytes=p.stat().st_size, sha256=digest(p)) for p in (calls, seqs)},
        validation=dict(complete_call_records=True, complete_sample_frames=True,
                        sample_count_matches_call_flags=True, sampled_calls=sample_count,
                        sampled_kinds=dict(sample_kinds), cigar_consumption_beyond_input=invalid_calls),
        total=present(total, total[1]),
        by_kind={k: present(s, total[1]) for k, s in zip(KINDS, by_kind) if s[0]},
        by_kind_and_max_input_length={k: {label: present(s, total[1]) for label, s in zip(LABELS, bins) if s[0]}
                                    for k, bins in zip(KINDS, by_size)},
        by_kind_and_flags={f'{KINDS[k]} flags={flags}': present(s, total[1]) for (k, flags), s in sorted(by_flags.items())},
        longest_calls=[dict(record_index=i, kind=KINDS[r[0]], milliseconds=ns / 1e6,
                            query_length=r[7], target_length=r[8], outer_bandwidth=r[4],
                            dp_flags=r[3], flags_hi=r[1]) for ns, i, r in sorted(slow, reverse=True)],
        caveats=[
            'Durations are summed per-helper wall times across workers, not elapsed runtime or hardware CPU samples.',
            'Instrumentation wraps alignment helpers. A helper can contain several probes, narrow attempts, subanchor segments, or rescue kernels.',
            'Recorded bandwidth is the outer caller parameter, not necessarily the actual inner SIMD bandwidth.',
            'The census does not wrap every path in the extension implementation; it is not exhaustive accounting of DP kernel invocations.',
            'Instrumented elapsed runtime is not a production speed benchmark.',
            'Sequence snapshots preserve helper inputs and outer parameters; exact inner-kernel replay requires the actual probe and retry arguments.',
            'Sample global sequence IDs are not used as unique identifiers because their generation uses a racing atomic load.'
        ])
    a.output.write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps({k: report[k] for k in ('validation', 'total', 'by_kind', 'by_kind_and_max_input_length')}, indent=2))

if __name__ == '__main__':
    main()
