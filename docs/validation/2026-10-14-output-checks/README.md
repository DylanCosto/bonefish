# Output correctness, determinism, input robustness and limits (14 Oct 2026)

Checks that a downstream tool or a clinical pipeline relies on, beyond accuracy: the files pass
strict validators, the output never depends on how bonefish was run, bad input fails loudly, and
large or unusual references work.

## 1. File-format conformance

Three bonefish outputs and minimap2 2.31 for comparison (`evidence/conformance/conformance.sh`):

| | input | output |
|---|---|---|
| A | genome-wide HG001 sample, FASTQ (16,177 reads) | sorted BAM, `--cs` |
| B | dorado BAM with `MM`/`ML`/`MN`/`RG` and 13 other dorado tags (20,000 HG002 chr20 reads) | sorted BAM, `-y --MD` |
| C | the same | sorted BAM, `-y -Y --cs` (soft-clipped supplementary records) |
| M | the same, minimap2 `-ax lr:hq -y --MD` + samtools sort | for comparison |

Each was checked with samtools quickcheck; Picard ValidateSamFile 3.5.0 with all checks, and again
with the reference on primary and supplementary records (which verifies NM and MD; Picard cannot take
a reference with the secondary records both aligners write without a sequence); samtools calmd, which
recomputes NM and MD from the reference; `check_sam.py`, which checks every record and read; and the
index (idxstats totals, 200 random region queries against a full scan).

| check | A | B | C | minimap2 |
|---|---|---|---|---|
| Picard, with reference (primary + supplementary) | missing read group only¹ | **no errors** | **no errors** | 5 warnings (adjacent I and D in a CIGAR) |
| Picard, all records | missing read group only¹; 27 warnings² | 53 warnings² | 53 warnings² | 7 CIGAR + 52 warnings² |
| NM / MD differences from the reference (calmd) | 0 / 0 | 0 / 0 | 0 / 0 | 0 / 0 |
| CIGAR length = sequence length | 16,948 records, 0 errors | 29,004, 0 | 29,004, 0 | 27,524, 0 |
| one primary record per read | 16,177, 0 | 20,000, 0 | 20,000, 0 | 20,000, 0 |
| SA tags match the records they point to | 2,010, 0 | 13,754, 0 | 13,754, 0 | 11,858, 0 |
| NM equals the edit distance in `cs`; `cs` matches the CIGAR | 20,490, 0 | | 44,032, 0 | |
| RG tags resolve to an `@RG` line | | 44,032, 0 | 44,032, 0 | 55,312, 0 |
| `MM`/`ML` valid for the record's sequence, `MN` consistent | | 29,004, 0 | 29,004, 0 | 27,524, 0 |
| index: idxstats total; region queries | exact; 200/200 | exact; 200/200 | exact; 200/200 | exact; 200/200 |

¹ FASTQ input has no read group unless `-R` is given (GATK-style tools require one; minimap2 is the
same). ² `QUALITY_NOT_STORED`: secondary records without a sequence, as minimap2 writes them.

- bonefish copies dorado's `@RG` line into the header and the `RG` tag onto every record; minimap2
  `-y` copies the tag but not the header line (it was added by hand for the comparison).
- Base-modification tags stay valid: on hard-clipped supplementary records `MN` differs from the
  sequence length, so modkit and samtools know the tags no longer apply (with `-Y` the full sequence
  is kept and the tags apply).
- With an aligned dorado BAM as input, bonefish drops the old alignment tags. Converting the same
  BAM with `samtools fastq -T '*'` and aligning with minimap2 `-y` gives records with two `NM` tags.

## 2. Determinism

The same HG001 sample (`evidence/determinism/determinism.sh`), SAM output hashed without `@PG`:

| varied | runs | distinct outputs |
|---|---|---:|
| threads | 1, 2, 4, 8, 16 | 1 |
| repeated runs | 3 at 16 threads | 1 |
| batch size `-K` | 5M, 50M, 500M | 1 |
| input form | FASTQ, gzip, bgzip, unaligned BAM, two files, stdin FASTQ, stdin BAM, pipe | 1 |
| index | raw, `.idx`, built from the FASTA | 1 |
| SIMD code path | AVX2, SSE4.1, scalar, automatic | 1 |

- Sorted BAM at 4 and 16 threads: identical records. Sorted BAM against SAM output: the same 21,089
  records (floats compared as numbers: SAM text has `de:f:0.1750` where samtools prints `0.175`).
- The genome-wide sample (161,762 reads) reproduces the published MD5 `deed7949d35f15bd8eed935be47a1d81`
  (v0.1.1 on x86, and on Apple M4 Pro), with this release's code and with its new default batch size.
- Thread counts 1-16 on the genome-wide sample: same MD5 ([scaling](../2026-10-14-scaling-memory/README.md)).

## 3. Input robustness

- **A bug found and fixed:** reads given through a pipe (`bonefish ref.fa <(zcat reads.fq.gz)`) failed
  with "not FASTA, FASTQ or BAM", and a dorado BAM streamed on stdin lost its `@RG` header lines. The
  pass that reads `@RG` lines opened the input separately; for a pipe it consumed the first bytes. Now
  the bytes read are replayed. Regression test `stdin_and_pipe_input`.
- **Malformed input** (`bonefish-core/tests/parser_stress.rs`): 3 million random and mutated inputs
  (truncations, bit flips, inserted line breaks, corrupted BGZF blocks and BAM fields, extreme length
  fields) through the FASTA, FASTQ, gzip, BGZF and BAM readers. It found 7 ways corrupt BAM/BGZF input
  could crash bonefish, and one silent failure: a corrupt bgzipped FASTQ was read as a shorter file
  with no error. All now give an error message and a non-zero exit; each case has a regression test.
  A decompression thread that stops early is now an error rather than the end of the input. After
  the fixes: no crash or hang in 3 million cases (and 1.5 million with overflow checks on).

## 4. Large references

Synthetic 4.6 Gbp references built from GRCh38 (`evidence/large_reference/large_ref.py`); the HG001
sample mapped against each and against GRCh38 with the same repeat threshold, and compared after
translating coordinates back:

| reference | what it tests | result |
|---|---|---|
| chr1 preceded by 1.5 Gbp of N, nothing else changed (195 contigs) | positions past 2^32 in the concatenated reference (chr19 onwards: 2,509 records); a 1.75 Gbp contig | all 16,132 reads not on chr1 identical; 45 reads on chr1 differ in small ways³ |
| 150,000 random 10 kb contigs, then GRCh38 with chr1-3 fused (689 Mbp) | an index of 150,193 contigs; positions past 2^32 | index in one part (87 s, 15 GB); no read placed on a decoy; 97.6% of reads identical⁴ |
| both: sorted BAM | contigs longer than BAI can address (2^29) | `.csi` index written; region queries on the last 89 Mbp of the long contig return every record (493/493, 491/491) |

- **A gap found and fixed:** sorted BAM output always wrote a BAI index, which cannot address positions
  past 512 Mbp (many plant chromosomes are longer). It now writes a CSI index when any reference is
  longer than 2^29 bases; human references still get the identical `.bai`.
- ³ Alignment heuristics inherited from minimap2 depend slightly on absolute positions: shifting chr1
  by just 1 kb changes 47 reads with bonefish and 50 with minimap2 (`shift_test.txt`). It does not
  grow with the offset (a 900 Mb shift leaves the read with the largest difference unchanged), so it is
  not an overflow.
- ⁴ Random decoys contain chance copies of about 0.5% of real 19-mers and move every chromosome's
  contig number, which minimap2-style tie-breaking uses, so some alignments legitimately change.

## Files

- `evidence/conformance/`: `summary.md` (all results), Picard reports, `check_sam.py` tables,
  `calmd_counts.tsv`, the dorado `@RG` line, `conformance.sh`, `check_sam.py`.
- `evidence/determinism/`: `results.tsv`, `determinism.sh`.
- `evidence/large_reference/`: per reference `summary.txt` and GNU time of the index build and
  mapping; `shift_test.txt`; `REFERENCES.sha256` (the generated references, which `large_ref.py`
  recreates exactly); the scripts.
- `SHA256SUMS` covers `evidence/`. Local paths are replaced by `${REPO}`, `${DATA}` and `${HOME}`.
