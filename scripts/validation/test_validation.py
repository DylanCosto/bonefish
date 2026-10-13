"""Small preservation and corrupt-input checks; no external dataset download."""
import gzip
import json
from pathlib import Path
import struct
import subprocess
import sys
import tempfile
import unittest

import analyze_census
from benchmark_sam import fastq_inventory, timing_metrics, validate_sam
from run_hg002_wgs import reference_contigs


class SamValidationTests(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.directory.cleanup)
        self.root = Path(self.directory.name)
        self.fastq = self.root / 'reads.fq.gz'
        with gzip.open(self.fastq, 'wb') as f:
            f.write(b'@read1\nACGTA\n+\n!I#$%\n@read2\nTTN\n+\nIII\n')
        self.expected = fastq_inventory([self.fastq])
        self.sam = self.root / 'out.sam'
        self.header = b'@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:100\n'
        self.reverse = b'read1\t16\tchr1\t1\t60\t5M\t*\t0\t0\tTACGT\t%$#I!\n'
        self.unmapped = b'read2\t4\t*\t0\t0\t*\t*\t0\t0\tTTN\tIII\n'
        self.secondary = b'read1\t256\tchr1\t10\t0\t5M\t*\t0\t0\t*\t*\n'

    def write(self, *records):
        self.sam.write_bytes(self.header + b''.join(records))

    def test_reverse_quality_unmapped_inventory_and_pg_normalization(self):
        self.write(self.reverse, self.unmapped, self.secondary)
        result = validate_sam(self.sam, self.expected)
        self.assertEqual((result['primary'], result['mapped_primary'], result['records']), (2, 1, 3))
        self.sam.write_bytes(b'@PG\tID:changed-command-path\n' + self.sam.read_bytes())
        self.assertEqual(validate_sam(self.sam, self.expected)['sam_sha256_without_pg'], result['sam_sha256_without_pg'])

    def test_missing_primary_is_rejected(self):
        self.write(self.reverse, self.secondary)
        with self.assertRaisesRegex(ValueError, 'Missing primary'):
            validate_sam(self.sam, self.expected)

    def test_duplicate_primary_is_rejected(self):
        self.write(self.reverse, self.unmapped, self.unmapped)
        with self.assertRaisesRegex(ValueError, 'Duplicate primary'):
            validate_sam(self.sam, self.expected)

    def test_changed_sequence_or_quality_is_rejected(self):
        for record in (self.reverse.replace(b'TACGT', b'TACGA'), self.reverse.replace(b'%$#I!', b'%$#II')):
            self.write(record, self.unmapped)
            with self.assertRaisesRegex(ValueError, 'differs from input'):
                validate_sam(self.sam, self.expected)

    def test_unexpected_secondary_read_is_rejected(self):
        self.write(self.reverse, self.unmapped, self.secondary.replace(b'read1', b'unknown'))
        with self.assertRaisesRegex(ValueError, 'unexpected read ID'):
            validate_sam(self.sam, self.expected)

    def test_duplicate_and_truncated_fastq_are_rejected(self):
        with self.assertRaisesRegex(ValueError, 'Duplicate FASTQ'):
            fastq_inventory([self.fastq, self.fastq])
        broken = self.root / 'truncated.fq'
        broken.write_bytes(b'@read1\nACGT\n+\n')
        with self.assertRaisesRegex(ValueError, 'FASTQ framing'):
            fastq_inventory([broken])

    def test_bsd_and_gnu_rss_units(self):
        self.assertEqual(timing_metrics('  12345 maximum resident set size\n')['peak_rss_bytes'], 12345)
        self.assertEqual(timing_metrics('Maximum resident set size (kbytes): 12345\n')['peak_rss_bytes'], 12345 * 1024)
        self.assertIsNone(timing_metrics('')['peak_rss_bytes'])

    def test_reference_lengths_from_plain_and_ordinary_gzip_fasta(self):
        plain = self.root / 'reference.fa'
        zipped = self.root / 'reference.fa.gz'
        plain.write_bytes(b'>chr1 description\nACGT\nNN\n>chr2\nACG\n')
        with gzip.open(zipped, 'wb') as f:
            f.write(plain.read_bytes())
        self.assertEqual(reference_contigs(plain), [('chr1', 6), ('chr2', 3)])
        self.assertEqual(reference_contigs(zipped), reference_contigs(plain))

    def test_duplicate_or_empty_reference_contigs_are_rejected(self):
        reference = self.root / 'reference.fa'
        for content in (b'>chr1\nA\n>chr1\nC\n', b'>chr1\n'):
            reference.write_bytes(content)
            with self.assertRaises(ValueError):
                reference_contigs(reference)


class CensusValidationTests(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.directory.cleanup)
        self.root = Path(self.directory.name)
        self.prefix = self.root / 'profile'
        rows = []
        for i in range(200):
            row = [0] * 19
            row[1], row[7], row[8], row[18] = (4 if i == 0 else 0), 3, 4, 100
            rows.append(analyze_census.REC.pack(*row))
        self.calls = self.root / 'profile.calls.bin'
        self.seqs = self.root / 'profile.seqs.bin'
        self.calls.write_bytes(b''.join(rows))
        self.seqs.write_bytes(struct.pack('<IQ', 0, 1) + rows[0] + bytes(29)
                              + struct.pack('<II', 3, 4) + bytes([0, 1, 4, 0, 1, 2, 3]))

    def run_decoder(self):
        return subprocess.run([sys.executable, '-B', str(Path(analyze_census.__file__)), str(self.prefix),
                               str(self.root / 'summary.json')], capture_output=True)

    def test_valid_fixed_and_variable_frames(self):
        result = self.run_decoder()
        self.assertEqual(result.returncode, 0, result.stderr.decode())
        summary = json.loads((self.root / 'summary.json').read_text())
        self.assertEqual(summary['total']['calls'], 200)
        self.assertEqual(summary['validation']['sampled_calls'], 1)

    def test_truncated_call_record_is_rejected(self):
        self.calls.write_bytes(self.calls.read_bytes()[:-1])
        self.assertNotEqual(self.run_decoder().returncode, 0)

    def test_truncated_sample_frame_is_rejected(self):
        self.seqs.write_bytes(self.seqs.read_bytes()[:-1])
        self.assertNotEqual(self.run_decoder().returncode, 0)

    def test_invalid_nt4_sample_is_rejected(self):
        self.seqs.write_bytes(self.seqs.read_bytes()[:-1] + b'\x05')
        self.assertNotEqual(self.run_decoder().returncode, 0)


if __name__ == '__main__':
    unittest.main()
