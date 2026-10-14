#!/usr/bin/env bash
# Output-format conformance of bonefish BAM/SAM:
#   A. genome-wide HG001 sample (16,177 reads) -> sorted BAM with --cs
#   B. dorado BAM with MM/ML/MN/RG tags (HG002 chr20, PAW70337) with -y --MD -> sorted BAM
#   C. as B with -Y --cs (soft-clipped supplementary records, so MM/ML stay valid on them)
#   (like minimap2, bonefish writes cs or MD, not both)
# Checks: samtools quickcheck; Picard ValidateSamFile (all checks; again with the reference on
# primary and supplementary records, so NM/MD are verified: Picard fails on secondary records
# without SEQ, which minimap2 and bonefish both write, when given a reference); samtools calmd (recomputes NM and MD from the reference); bench/check_sam.py
# (CIGAR/SEQ, SA cross-references, NM vs cs, RG, MAPQ, MM/ML/MN); the .bai index (idxstats total
# and 200 random region queries vs a full scan). minimap2 2.31 output for B is checked the same way.
# Usage: bench/conformance.sh OUTDIR     Env: BONEFISH
set -uo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
export PATH="${MAMBA_ROOT_PREFIX:-$HOME/micromamba}/envs/vareval/bin:$PATH"
PICARD="${MAMBA_ROOT_PREFIX:-$HOME/micromamba}/envs/picard/bin/picard"
MM="${MAMBA_ROOT_PREFIX:-$HOME/micromamba}/envs/ontaln/bin/minimap2"
B="${BONEFISH:-$ROOT/bonefish/target/release/bonefish}"
OUT="${1:?outdir}"; W="$OUT/work"; mkdir -p "$W"
REF="$ROOT/data/ref/GRCh38_no_alt.fa"; RAW="$ROOT/data/index/bonefish.lrhq.raw"; MMI="$ROOT/data/index/mm2.lrhq.mmi"
DICT="${REF%.fa}.dict"
[[ -s "$DICT" ]] || samtools dict -o "$DICT" "$REF"
S="$OUT/summary.md"; : > "$S"

# dorado input: the first 20,000 primary records from chr20:30 Mb of the dorado BAM (a sparse
# local copy holding only the chr20 blocks, so it is read by region, single-threaded:
# samtools view -@ hangs on it)
DOR="$W/dorado20k.bam"
if [[ ! -s "$DOR" ]]; then
  samtools view -h -F 0x900 "$ROOT/data/reads/PAW70337.sparse.bam" chr20:30000000-64444167 2>/dev/null \
    | awk '/^@/ || ++n <= 20000' | samtools view -b -o "$DOR" -
fi
echo "dorado input: $(samtools view -c "$DOR") records" | tee -a "$S"
samtools view -H "$DOR" | grep '^@RG' | cut -c1-150 > "$OUT/input_rg.txt"

cat "$RAW" > /dev/null
[[ -s "$W/A.bam" ]] || "$B" -t 16 --cs -o "$W/A.bam" "$RAW" "$ROOT/data/reads/hg001_wg_sub10.fq" 2> "$W/A.log"
[[ -s "$W/B.bam" ]] || "$B" -t 16 -y --MD -o "$W/B.bam" "$RAW" "$DOR" 2> "$W/B.log"
[[ -s "$W/C.bam" ]] || "$B" -t 16 -y -Y --cs -o "$W/C.bam" "$RAW" "$DOR" 2> "$W/C.log"
if [[ ! -s "$W/M.bam" ]]; then
  # dorado's own tags only (the input's alignment tags would be duplicated by minimap2 -y)
  samtools fastq -T MM,ML,MN,RG,qs,du,ns,ts,mx,ch,st,rn,fn,sm,sd,sv,dx,pi,sp,pt,bh "$DOR" 2>/dev/null > "$W/dor.fq"
  "$MM" -ax lr:hq -y --MD -t 16 "$MMI" "$W/dor.fq" 2> "$W/M.log" | samtools sort -@ 8 -o "$W/M.bam" --write-index -
  # minimap2 does not copy the input @RG line: add it so RG tags resolve
  samtools view -H "$W/M.bam" > "$W/M.hdr"; grep '^@RG' <(samtools view -H "$DOR") >> "$W/M.hdr"
  samtools reheader "$W/M.hdr" "$W/M.bam" > "$W/M2.bam" && mv "$W/M2.bam" "$W/M.bam" && rm -f "$W/M.bam.csi" && samtools index "$W/M.bam"
  rm -f "$W/dor.fq"
fi

for x in A B C M; do
  bam="$W/$x.bam"
  echo "### $x: $(samtools view -c "$bam") records" | tee -a "$S"
  q=$(samtools quickcheck -vvv "$bam" 2>&1 | tail -1); echo "- samtools quickcheck: ${q:-ok}" | tee -a "$S"
  "$PICARD" ValidateSamFile -I "$bam" -M SUMMARY > "$W/$x.picard.txt" 2> "$W/$x.picard.log"
  echo "- Picard ValidateSamFile (exit $?): $(grep -E 'ERROR|WARNING|No errors' "$W/$x.picard.txt" | tr '\n' ' ' | tr -s ' ')" | tee -a "$S"
  samtools view -b -F 0x100 -o "$W/$x.nosec.bam" "$bam"
  "$PICARD" ValidateSamFile -I "$W/$x.nosec.bam" -R "$REF" -M SUMMARY > "$W/$x.picard_ref.txt" 2> "$W/$x.picard_ref.log"
  echo "- Picard ValidateSamFile with reference, primary + supplementary (exit $?): $(grep -E 'ERROR|WARNING|No errors' "$W/$x.picard_ref.txt" | tr '\n' ' ' | tr -s ' ')" | tee -a "$S"
  rm -f "$W/$x.nosec.bam"
  samtools calmd -@ 8 "$bam" "$REF" 2> "$W/$x.calmd.err" > /dev/null
  echo "- samtools calmd: $(grep -c 'different NM' "$W/$x.calmd.err") NM and $(grep -c 'different MD' "$W/$x.calmd.err") MD differences from the reference" | tee -a "$S"
  samtools view -h "$bam" | python3 "$ROOT/bench/check_sam.py" "$x" > "$W/$x.check.md"; echo "- check_sam.py exit $?" | tee -a "$S"
  cat "$W/$x.check.md" >> "$S"
  # index: idxstats total = record count; random regions = full-scan counts
  tot=$(samtools idxstats "$bam" | awk '{s+=$3+$4} END {print s}')
  echo "- idxstats total $tot vs $(samtools view -c "$bam") records" | tee -a "$S"
  samtools view "$bam" | awk '!and($2,4){print $3"\t"$4}' > "$W/$x.pos"
  python3 - "$bam" "$W/$x.pos" <<'PY' | tee -a "$S"
import random, subprocess, sys, collections
bam, posf = sys.argv[1], sys.argv[2]
starts = collections.defaultdict(list)
for line in open(posf):
    c, p = line.split(); starts[c].append(int(p))
contigs = [c for c in starts if len(starts[c]) > 0]
rnd = random.Random(7); bad = 0
for _ in range(200):
    c = rnd.choice(contigs); p = rnd.choice(starts[c]); a, b = p, p + rnd.choice([1, 1000, 100000])
    # records starting inside [a, b] (region queries also return overlapping records; compare starts)
    out = subprocess.run(["samtools", "view", bam, f"{c}:{a}-{b}"], capture_output=True, text=True).stdout
    got = sum(1 for l in out.splitlines() if a <= int(l.split("\t")[3]) <= b and not int(l.split("\t")[1]) & 4)
    want = sum(1 for s in starts[c] if a <= s <= b)
    bad += got != want
print(f"- index region queries: {200 - bad}/200 agree with a full scan")
PY
  rm -f "$W/$x.pos"
done
