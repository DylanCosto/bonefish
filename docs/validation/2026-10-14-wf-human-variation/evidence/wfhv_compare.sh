#!/usr/bin/env bash
# ONT's wf-human-variation (v2.8.1) on the same HG002 chr20 reads aligned by minimap2 and by bonefish
# (results/modkit/*.bam, from bench/modkit_compare.sh): small variants (Clair3), SVs (Sniffles2),
# methylation (modkit), STRs (Straglr) and CNVs, restricted to chr20.
set -uo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
W="$ROOT/results/wfhv"; M="$ROOT/results/modkit"
# Work and output directories on a path without spaces: Clair3 passes unquoted paths to gzip, and the data
# volume is mounted at "<data volume mount point>". Outputs are copied to $W afterwards; work dirs are deleted.
L="$HOME/wfhv_work"; mkdir -p "$L"
NF=$HOME/micromamba/envs/nf/bin/nextflow
export NXF_HOME="$L/.nextflow" NXF_ANSI_LOG=false NXF_SYNTAX_PARSER=v1  # v2.8.1 predates the strict parser of Nextflow 26
log() { echo "[$(date +%T)] $*" | tee -a "$W/run.log"; }
cd "$W"
for a in minimap2 bonefish; do
  [[ -f $W/$a/.done ]] && continue
  log "start $a"
  /usr/bin/time -v -o $W/$a.time $NF run epi2me-labs/wf-human-variation -r v2.8.1 -c $W/local.config -w $L/work_$a \
    --bam $M/$a.bam --ref $ROOT/data/ref/GRCh38_no_alt.fa --bed $W/chr20.bed \
    --snp --sv --mod --str --cnv --sex XY --override_basecaller_cfg dna_r10.4.1_e8.2_400bps_sup@v5.0.0 \
    --sample_name HG002_$a --threads 16 --out_dir $L/out_$a > $W/$a.nf.log 2>&1
  rc=$?; log "end $a rc=$rc"
  if [[ $rc == 0 ]]; then
    mkdir -p $W/$a && cp -rL $L/out_$a/. $W/$a/ && touch $W/$a/.done && rm -rf $L/work_$a $L/out_$a
  fi
done
log ALLDONE
