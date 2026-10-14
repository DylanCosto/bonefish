#!/usr/bin/env bash
# Speed and mapping statistics of ONT aligners on the same reads, reference and machine.
#   1. index builds (timed);
#   2. SAM output on a 10% genome-wide subset (all tools) and on the 1x HG001 sample (fast tools);
#   3. reads -> sorted, indexed BAM on the 1x sample (bonefish, dorado aligner, minimap2 + samtools).
# Each run: GNU time -v (wall, CPU, peak RSS); SAM runs keep columns 1-6,12- for bench/bam_summary.py.
# Tools: minimap2 (ontaln env); winnowmap, meryl, ngmlr, lra, mm2plus (alnbench env); mm2-fast built
# from source; dorado 2.1.2 (ONT release). Usage: bench/aligner_bench.sh [OUT_DIR]
set -uo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
OUT="${1:-$ROOT/results/alnbench}"; mkdir -p "$OUT/runs" "$OUT/stats"
T=16
REF="$ROOT/data/ref/GRCh38_no_alt.fa"
IDX="$ROOT/data/index/bench"; mkdir -p "$IDX"
READS="$ROOT/data/reads/hg001_wg_1x.fq.gz"
SUB="$ROOT/data/reads/hg001_wg_sub10.fq"
E=$HOME/micromamba/envs/alnbench/bin
MM2=$HOME/micromamba/envs/ontaln/bin/minimap2
SAMTOOLS=$HOME/micromamba/envs/vareval/bin/samtools
TOOLS="${DATA}/tools"
MM2FAST="$TOOLS/mm2-fast_src/mm2-fast/minimap2"
DORADO="$TOOLS/dorado-2.1.2-linux-x64/bin/dorado"
BF="$ROOT/bonefish/target/release/bonefish"
export LD_LIBRARY_PATH="$HOME/micromamba/envs/alnbench/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
STAGE="$OUT/stage"; echo setup > "$STAGE"
MON_CSV="$OUT/resources.csv" MON_STAGE="$STAGE" "$ROOT/bench/monitor.sh" 10 60 > "$OUT/monitor.log" 2>&1 &
MON=$!
trap 'kill $MON 2>/dev/null' EXIT
log() { echo "[$(date +%T)] $*" | tee -a "$OUT/driver.log"; }

# timed NAME CMD...: run CMD (a bash string) with GNU time -v on the whole pipeline's first command
timed() {
  local name=$1; shift
  echo "$name" > "$STAGE"; log "start $name"
  bash -c "$*" 2> "$OUT/runs/$name.err"; local rc=$?
  log "end $name rc=$rc $(grep -E 'Elapsed|Maximum resident' "$OUT/runs/$name.time" 2>/dev/null | sed 's/^\s*//' | tr '\n' ' ')"
}
TV="/usr/bin/time -v -o"

{ echo "bonefish $($BF --version)"; echo "minimap2 $($MM2 --version)"; echo "mm2plus $($E/mm2plus.avx2 --version 2>/dev/null)";
  echo "mm2-fast $($MM2FAST --version 2>/dev/null | tail -1)"; echo "dorado $("$DORADO" --version 2>&1 | tail -1)";
  echo "winnowmap $($E/winnowmap --version)"; echo "ngmlr $($E/ngmlr --version 2>&1 | head -1)"; echo "lra $($E/lra 2>&1 | grep Version)"; } > "$OUT/versions.txt"

# --- subset: every 10th read of the 1x sample ---
if [[ ! -s "$SUB" ]]; then
  log "making subset"
  zcat "$READS" | awk 'NR%40>=1 && NR%40<=4' > "$SUB"
fi

# --- 1. indexes ---
[[ -s $IDX/bonefish.raw ]] || timed idx_bonefish "$TV $OUT/runs/idx_bonefish.time $BF -t $T -d $IDX/bonefish.raw --idx-raw $REF"
[[ -s $IDX/mm2.lrhq.mmi ]] || timed idx_minimap2 "$TV $OUT/runs/idx_minimap2.time $MM2 -t $T -x lr:hq -d $IDX/mm2.lrhq.mmi $REF"
[[ -s $IDX/mm2plus.lrhq.mmi ]] || timed idx_mm2plus "$TV $OUT/runs/idx_mm2plus.time $E/mm2plus.avx2 -t $T -x lr:hq -d $IDX/mm2plus.lrhq.mmi $REF"
[[ -s $IDX/mm2fast.mapont.mmi ]] || timed idx_mm2fast "$TV $OUT/runs/idx_mm2fast.time '$MM2FAST' -t $T -x map-ont -d $IDX/mm2fast.mapont.mmi $REF"
if [[ ! -s $IDX/winnowmap.rep15.txt ]]; then
  timed idx_meryl_count "$TV $OUT/runs/idx_meryl_count.time $E/meryl count k=15 threads=$T output $IDX/merylDB_k15 $REF"
  timed idx_meryl_print "$TV $OUT/runs/idx_meryl_print.time $E/meryl print greater-than distinct=0.9998 $IDX/merylDB_k15 > $IDX/winnowmap.rep15.txt"
fi
# Winnowmap 2.03 maps nothing from a prebuilt index (checked on test data), so it indexes the FASTA on
# every run; this run on 10 reads measures that index build, to subtract from its mapping times.
head -40 "$SUB" > "$IDX/tiny.fq"
[[ -s $OUT/runs/idx_winnowmap.time ]] || timed idx_winnowmap "$TV $OUT/runs/idx_winnowmap.time $E/winnowmap -t $T -W $IDX/winnowmap.rep15.txt -ax map-ont $REF $IDX/tiny.fq > /dev/null"
if [[ ! -s $REF-ht-13-2.2.ngm ]]; then
  timed idx_ngmlr "$TV $OUT/runs/idx_ngmlr.time $E/ngmlr -x ont -t $T -r $REF -q $IDX/tiny.fq -o /dev/null"
fi
[[ -s $REF.gli ]] || timed idx_lra "$TV $OUT/runs/idx_lra.time $E/lra index -ONT $REF"

# --- 2. SAM output ---
# sam NAME READS INDEXFILE CMD: warm the index in the page cache, then run CMD > SAM columns
sam() {
  local name=$1 reads=$2 warm=$3; shift 3
  [[ -s "$OUT/stats/$name/bam_summary.tsv" ]] && return
  [[ -n "$warm" ]] && cat "$warm" > /dev/null
  timed "$name" "$TV $OUT/runs/$name.time $* | cut -f1-6,12- > $OUT/runs/$name.cols"
  mkdir -p "$OUT/stats/$name"
  grep -v '^@' "$OUT/runs/$name.cols" | python3 "$ROOT/bench/bam_summary.py" "$OUT/stats/$name"
  rm -f "$OUT/runs/$name.cols"
}
for set in sub10 1x; do
  R=$SUB; [[ $set == 1x ]] && R=$READS
  sam ${set}_bonefish "$R" $IDX/bonefish.raw "$BF -t $T -a $IDX/bonefish.raw $R"
  sam ${set}_minimap2 "$R" $IDX/mm2.lrhq.mmi "$MM2 -t $T -ax lr:hq $IDX/mm2.lrhq.mmi $R"
  sam ${set}_mm2plus "$R" $IDX/mm2plus.lrhq.mmi "$E/mm2plus.avx2 -t $T -ax lr:hq $IDX/mm2plus.lrhq.mmi $R"
  sam ${set}_mm2fast "$R" $IDX/mm2fast.mapont.mmi "'$MM2FAST' -t $T -ax map-ont $IDX/mm2fast.mapont.mmi $R"
  sam ${set}_dorado "$R" $IDX/mm2.lrhq.mmi "'$DORADO' aligner -t $T --emit-sam $IDX/mm2.lrhq.mmi $R"
  sam ${set}_winnowmap "$R" "" "$E/winnowmap -t $T -W $IDX/winnowmap.rep15.txt -ax map-ont $REF $R"
  if [[ $set == sub10 ]]; then
    sam ${set}_ngmlr "$R" "" "$E/ngmlr -x ont -t $T -r $REF -q $R"
    sam ${set}_lra "$R" "" "$E/lra align -ONT -t $T $REF $R -p s"
  fi
done

# --- 3. reads -> sorted, indexed BAM (1x) ---
B="$OUT/bam"; mkdir -p "$B"
[[ -s $B/bonefish.bam.bai ]] || { cat $IDX/bonefish.raw > /dev/null; timed bam_bonefish "$TV $OUT/runs/bam_bonefish.time $BF -t $T -o $B/bonefish.bam $IDX/bonefish.raw $READS"; }
[[ -d $B/dorado ]] || { cat $IDX/mm2.lrhq.mmi > /dev/null; timed bam_dorado "$TV $OUT/runs/bam_dorado.time '$DORADO' aligner -t $T -o $B/dorado $IDX/mm2.lrhq.mmi $READS"; }
[[ -s $B/minimap2.bam.csi || -s $B/minimap2.bam.bai ]] || { cat $IDX/mm2.lrhq.mmi > /dev/null; timed bam_minimap2 "$TV $OUT/runs/bam_minimap2.time bash -c \"$MM2 -t $T -ax lr:hq $IDX/mm2.lrhq.mmi $READS | $SAMTOOLS sort -@$T -m 384M --write-index -o $B/minimap2.bam -\""; }
echo done > "$STAGE"; log ALLDONE
