#!/usr/bin/env bash
# Limits test on 4.6 Gbp synthetic references (bench/large_ref.py). The HG001 sample (16,177 reads)
# is mapped against the synthetic reference and against plain GRCh38 with the same fixed repeat
# threshold (-f), and the outputs are compared after translating coordinates back to GRCh38.
#   DECOY=pad (default): 1.5 Gbp of N before chr1's sequence, nothing else changed. Contig order is
#     unchanged (minimap2-style aligners break ties with it), so output must be identical for every
#     read not on chr1; chr1 is 1.75 Gbp (past the BAI limit of 2^29) and chr19 onwards lie past 2^32
#     in the concatenated reference.
#   DECOY=random: 150,000 random 10 kb contigs first, then GRCh38 with chr1-3 fused into one 689 Mbp
#     contig: many contigs, positions past 2^32, a contig past the BAI limit; the decoys add chance
#     seed matches and shift contig order, so some alignments legitimately differ.
# Also checks the sorted BAM and its index on the long contig.
# Usage: bench/large_ref.sh OUTDIR     Env: BONEFISH, DECOY
set -uo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
export PATH="${MAMBA_ROOT_PREFIX:-$HOME/micromamba}/envs/vareval/bin:$PATH"
B="${BONEFISH:-$ROOT/bonefish/target/release/bonefish}"
OUT="${1:?outdir}"; W="$ROOT/data/ref/extra/large"; mkdir -p "$OUT" "$W"
REF="$ROOT/data/ref/GRCh38_no_alt.fa"; RAW="$ROOT/data/index/bonefish.lrhq.raw"
Q="$ROOT/data/reads/hg001_wg_sub10.fq"
DECOY="${DECOY:-pad}"
case $DECOY in
  pad) SUF=P; NDEC=0; LONG=chr1 ;;
  random) SUF=""; NDEC=150000; LONG=chr1_2_3 ;;
  *) echo "DECOY must be pad or random"; exit 1 ;;
esac
BIG="$W/large$SUF.fa"; IX="$W/large$SUF.raw"
[[ -s "$BIG" ]] || python3 "$ROOT/bench/large_ref.py" "$REF" "$NDEC" "$BIG" "$DECOY"
[[ -s "$BIG.fai" ]] || samtools faidx "$BIG"
awk -v d="$DECOY" '{s+=$2; n++} END {printf "decoys: %s; reference: %d contigs, %.3f Gbp\n", d, n, s/1e9}' "$BIG.fai" | tee "$OUT/summary.txt"
LEN=$(awk -v c="$LONG" '$1==c{print $2}' "$BIG.fai"); echo "$LONG length $LEN" | tee -a "$OUT/summary.txt"

# repeat threshold of plain GRCh38 (printed when an index is built), fixed for both runs
MID=$("$B" -t 16 -a "$REF" "$ROOT/data/index/bench/tiny.fq" 2>&1 >/dev/null | sed -n 's/.*mid_occ: \([0-9]*\).*/\1/p' | head -1)
[[ -n "$MID" ]] || { echo "could not read the repeat threshold"; exit 1; }
echo "fixed repeat threshold -f $MID" | tee -a "$OUT/summary.txt"

[[ -s "$IX" ]] || /usr/bin/time -f "%e s wall, %M KB peak" -o "$OUT/index_build.time" "$B" -t 16 -d "$IX" --idx-raw "$BIG" 2> "$OUT/index_build.log"
[[ -s "$OUT/index_build.time" ]] && echo "index build: $(cat "$OUT/index_build.time"); $(grep -o 'Done in.*part)' "$OUT/index_build.log")" | tee -a "$OUT/summary.txt"
/usr/bin/time -f "%e s wall, %M KB peak" -o "$OUT/map_large.time" "$B" -a -t 16 -f "$MID" "$IX" "$Q" 2> "$OUT/map_large.log" > "$W/large.sam"
"$B" -a -t 16 -f "$MID" "$RAW" "$Q" 2> "$OUT/map_plain.log" > "$W/plain.sam"
echo "mapping vs large: $(cat "$OUT/map_large.time")" | tee -a "$OUT/summary.txt"

python3 - "$W/plain.sam" "$W/large.sam" "$DECOY" <<'PY' | tee -a "$OUT/summary.txt"
import collections, sys
plain, large, decoy = sys.argv[1:4]
lens = {}
for line in open(plain):
    if not line.startswith("@"): break
    if line.startswith("@SQ"):
        t = dict(x.split(":", 1) for x in line.rstrip().split("\t")[1:])
        lens[t["SN"]] = int(t["LN"])
b1, b2 = lens["chr1"], lens["chr1"] + lens["chr2"]
def tr(c, p):
    if decoy == "pad":
        return (c, p - 1_500_000_000) if c == "chr1" else (c, p)
    if c != "chr1_2_3": return c, p
    if p <= b1: return "chr1", p
    if p <= b2: return "chr2", p - b1
    return "chr3", p - b2
def norm(line):
    f = line.rstrip("\n").split("\t")
    if f[2] != "*":
        f[2], p = tr(f[2], int(f[3])); f[3] = str(p)
    for i in range(11, len(f)):
        if f[i].startswith("SA:Z:"):
            ents = []
            for e in f[i][5:].split(";"):
                if not e: continue
                x = e.split(","); x[0], p = tr(x[0], int(x[1])); x[1] = str(p); ents.append(",".join(x))
            f[i] = "SA:Z:" + ";".join(ents) + ";"
    return "\t".join(f)
a = [l.rstrip("\n") for l in open(plain) if not l.startswith("@")]
b = [norm(l) for l in open(large) if not l.startswith("@")]
ra, rb = collections.defaultdict(set), collections.defaultdict(set)
for l in a: ra[l.split("\t")[0]].add(l)
for l in b: rb[l.split("\t")[0]].add(l)
diff = [q for q in ra if ra[q] != rb.get(q)]
near = ("chr1",) if decoy == "pad" else ("chr1", "chr2", "chr3")
on = sum(1 for q in diff if any(x.split("\t")[2] in near or ("SA:Z:chr1" in x) for x in ra[q] | rb.get(q, set())))
past = sum(1 for l in b if l.split("\t")[2] in ("chr19", "chr20", "chr21", "chr22", "chrX", "chrY", "chrM"))
print(f"reads: {len(ra):,}; identical: {len(ra) - len(diff):,}; different: {len(diff):,}, of which "
      f"{on:,} involve {'/'.join(near)} and {len(diff) - on:,} do not")
print(f"records: plain {len(a):,}, large {len(b):,}; on decoys: {sum(1 for l in b if l.split(chr(9))[2].startswith('decoy_')):,}; "
      f"on contigs past 2^32 in the concatenated reference: {past:,}")
PY

# sorted BAM with the long contig, and its index
"$B" -t 16 -f "$MID" -o "$W/large.bam" "$IX" "$Q" 2> "$OUT/sorted.log"; echo "sorted BAM exit $?" | tee -a "$OUT/summary.txt"
ls "$W"/large.bam.* 2>/dev/null | sed 's|.*/||' | tr '\n' ' ' | sed 's/^/index files: /' | tee -a "$OUT/summary.txt"; echo | tee -a "$OUT/summary.txt"
samtools quickcheck -vvv "$W/large.bam" 2>&1 | tail -1 | sed 's|^.*/|quickcheck: |' | tee -a "$OUT/summary.txt"
FROM=$((LEN - 89445510))   # the last ~89 Mbp of the long contig, past 2^29
want=$(samtools view "$W/large.bam" | awk -v c="$LONG" -v f="$FROM" '$3==c && $4>f' | wc -l)
got=$(samtools view "$W/large.bam" "$LONG:$((FROM + 1))-$LEN" 2>&1 | awk -F'\t' -v f="$FROM" 'NF>10 && $4>f' | wc -l)
echo "records starting past $FROM on $LEN-bp $LONG: full scan $want, CSI region query $got" | tee -a "$OUT/summary.txt"
rm -f "$W/large.sam" "$W/plain.sam" "$W/large.bam" "$W/large.bam.csi" "$W/large.bam.bai"
