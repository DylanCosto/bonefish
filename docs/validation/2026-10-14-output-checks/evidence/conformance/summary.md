dorado input: 20000 records
### A: 21089 records
- samtools quickcheck: results/conformance/work/A.bam has good EOF block.
- Picard ValidateSamFile (exit 2): ERROR:MISSING_READ_GROUP	1 WARNING:QUALITY_NOT_STORED	27 WARNING:RECORD_MISSING_READ_GROUP	21089 
- Picard ValidateSamFile with reference, primary + supplementary (exit 2): ERROR:MISSING_READ_GROUP	1 WARNING:RECORD_MISSING_READ_GROUP	17547 
- samtools calmd: 0 NM and 0 MD differences from the reference
- check_sam.py exit 0
## A: 16,177 reads
| check | records examined | failures | examples |
|---|---:|---:|---|
| cigar_seq | 16,948 | 0 |  |
| hardclip | 1,370 | 0 |  |
| one_primary | 16,177 | 0 |  |
| sa_entries | 2,010 | 0 |  |
| nm_cs | 20,490 | 0 |  |
| rg | 0 | 0 |  |
| mapq | 21,089 | 0 |  |
| mm_ml | 0 | 0 |  |
- idxstats total 21089 vs 21089 records
- index region queries: 200/200 agree with a full scan
### B: 44033 records
- samtools quickcheck: results/conformance/work/B.bam has good EOF block.
- Picard ValidateSamFile (exit 1): WARNING:QUALITY_NOT_STORED	53 
- Picard ValidateSamFile with reference, primary + supplementary (exit 0): No errors found 
- samtools calmd: 0 NM and 0 MD differences from the reference
- check_sam.py exit 0
## B: 20,000 reads
| check | records examined | failures | examples |
|---|---:|---:|---|
| cigar_seq | 29,004 | 0 |  |
| hardclip | 9,005 | 0 |  |
| one_primary | 20,000 | 0 |  |
| sa_entries | 13,754 | 0 |  |
| nm_cs | 0 | 0 |  |
| rg | 44,032 | 0 |  |
| mapq | 44,033 | 0 |  |
| mm_ml | 29,004 | 0 |  |
- idxstats total 44033 vs 44033 records
- index region queries: 200/200 agree with a full scan
### C: 44033 records
- samtools quickcheck: results/conformance/work/C.bam has good EOF block.
- Picard ValidateSamFile (exit 1): WARNING:QUALITY_NOT_STORED	53 
- Picard ValidateSamFile with reference, primary + supplementary (exit 0): No errors found 
- samtools calmd: 0 NM and 0 MD differences from the reference
- check_sam.py exit 0
## C: 20,000 reads
| check | records examined | failures | examples |
|---|---:|---:|---|
| cigar_seq | 29,004 | 0 |  |
| hardclip | 0 | 0 |  |
| one_primary | 20,000 | 0 |  |
| sa_entries | 13,754 | 0 |  |
| nm_cs | 44,032 | 0 |  |
| rg | 44,032 | 0 |  |
| mapq | 44,033 | 0 |  |
| mm_ml | 29,004 | 0 |  |
- idxstats total 44033 vs 44033 records
- index region queries: 200/200 agree with a full scan
### M: 55312 records
- samtools quickcheck: results/conformance/work/M.bam has good EOF block.
- Picard ValidateSamFile (exit 1): WARNING:ADJACENT_INDEL_IN_CIGAR	7 WARNING:QUALITY_NOT_STORED	52 
- Picard ValidateSamFile with reference, primary + supplementary (exit 1): WARNING:ADJACENT_INDEL_IN_CIGAR	5 
- samtools calmd: 0 NM and 0 MD differences from the reference
- check_sam.py exit 0
## M: 20,000 reads
| check | records examined | failures | examples |
|---|---:|---:|---|
| cigar_seq | 27,524 | 0 |  |
| hardclip | 7,524 | 0 |  |
| one_primary | 20,000 | 0 |  |
| sa_entries | 11,858 | 0 |  |
| nm_cs | 0 | 0 |  |
| rg | 55,312 | 0 |  |
| mapq | 55,312 | 0 |  |
| mm_ml | 27,524 | 0 |  |
- idxstats total 55312 vs 55312 records
- index region queries: 200/200 agree with a full scan
