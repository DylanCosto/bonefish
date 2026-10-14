# HG001 genome-wide Mac parity run

Tested source: `c1b7eb074b49806a86bbd0e48c617957ac436618` (clean NEON DP).
Drive FASTQ: [hg001_wg_1x.fq.gz](https://drive.google.com/file/d/1L2O5CAHXbbypQxm52ZLbp6xAqCYYS7XH/view); compressed SHA256 `2edff9a3834d5d55963da198bad3717064e0ed5232a037d6783c84fbc8ffb288`.
The download matched Drive's SHA256, SHA1 and MD5; file size was 1,945,506,828 bytes.

The standard `bonefish -a -t 8 REFERENCE READS` run completed in
**88.172 s**, including fresh index and SAM writing.
FASTQ, reference and SAM were on the internal APFS SSD. This is one run, not a paired speed comparison.
Exact commands and input/binary hashes are in [the run receipt](evidence/candidate-1.json); [build context](evidence/BUILD_CONTEXT.json) records the compiler, profile and host.

Input: **161,762 reads / 2,713,864,134 bases**, read N50 28,334 bp.
Output: **212,292 records**, 155,784 mapped and
5,978 unmapped primaries. Every input primary, sequence and quality
was preserved after restoring orientation. Full CIGAR consumption, SEQ lengths,
reference headers and bounds passed for every alignment record.

| Hash normalization | Complete SHA256 |
|---|---|
| Raw SAM, including all headers | `a2cb5bf5f4e77c66d29e38f7f03482c25044f82be641c8ef451f498b1201887e` |
| SAM excluding only @PG | `8da232f13b35280925914bc039132dddfaf61dfbcbf9131f7fa9179073eec8f1` |
| Exact record stream excluding all headers | `3a3e23a17db4fea47e2d5abbfd713737eaa7aed09802603b6f6ab2f128e1c275` |

Expected x86 prefix supplied: **`deed7949d35f`**.
Matching normalization: **MD5 of SAM excluding only @PG**.
Complete Mac MD5: **`deed7949d35f15bd8eed935be47a1d81`**.
The supplied x86 prefix is reproduced; the expected algorithm is inferred to be MD5.
The full x86 MD5, exact invocation and reference receipt remain unprovided.
This is a successful reproduction of the supplied prefix; full-value x86/ARM confirmation
remains pending. SHA256 of the same @PG-filtered stream is independently recorded above.
The record-stream normalization excludes headers only; it preserves all records,
order, CIGARs, sequences, qualities and tags. No fields were removed to force a match.

Recompute the matching fingerprint on macOS from the retained SAM:

```bash
grep -v '^@PG' candidate-1.sam | md5
grep -v '^@PG' candidate-1.sam | shasum -a 256
```

Both commands retain every other header and every record field in original order.

The complete SAM is retained in internal scratch. Large FASTQ/reference/SAM files stay outside Git.
This sample alone does not reproduce the separate GIAB variant-calling benchmarks.

Additional MD5 normalizations are in [ADDITIONAL_HASHES.json](evidence/ADDITIONAL_HASHES.json); sorted record-stream hashes
are in [SORTED_RECORD_HASHES.json](evidence/SORTED_RECORD_HASHES.json). Both retain all fields and record multiplicity.


[Full results](evidence/PARITY_RESULTS.json) include all validation flags, counts and hashes.
[Download receipt](evidence/DOWNLOAD_STATUS.json) and [source metadata](evidence/DRIVE_SOURCE.json)
record the exact Drive ID and byte hashes. The [reference receipt](evidence/REFERENCE_RECEIPT.json)
records internal staging of the 195-contig GRCh38 no-alt analysis set, compressed SHA256
`3c8def6d325c5d1e934b2dd530c4d1709f027677a15859467071fddf3fff2026`.
The sample's identity is supplied by the user and Drive filename; no independent donor or chemistry
qualification was performed.

The Google Drive connector rejected this 1.95 GB file at its 256 MiB limit. The existing
rclone connection downloaded it in 63.0 s including hash verification. Acquisition and
verification are outside the aligner timer. The [download log](evidence/download.log),
[aligner log](evidence/candidate-1.log) and [driver receipt log](evidence/default-sam.driver.log)
are retained. No aligner source changes were made for this replay.

Local paths are replaced with documented placeholders. The [import manifest](evidence/IMPORT_MANIFEST.json)
retains original and packaged hashes; [SHA256SUMS](evidence/SHA256SUMS.txt) covers every evidence file.
The [analysis scripts](evidence/analysis/) are historical local snapshots with substituted paths.
The first script snapshot predates identification of the matching MD5; the additional MD5 results
and final report retain that finding without rewriting the original measurement script.

```bash
cd docs/validation/2026-10-13-hg001-macos-parity/evidence
shasum -a 256 -c SHA256SUMS.txt
```
