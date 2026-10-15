#!/usr/bin/env bash
# Runs the 14 Oct 2026 validation set in sequence (each step skips work already done):
# simulated reads, HAC / R9.4.1 / sup chr20 variant calling, phasing, CMRG + tandem repeats,
# determinism, large reference, then thread scaling last on an otherwise idle machine.
set -uo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
R="$ROOT/results"; G="$ROOT/data/gen"
log() { echo "[$(date '+%F %T')] $*"; }
# wait for the simulations and conformance runs started separately
while ps aux | grep -qE '[b]adread simulate|[s]im_wg\.sh|[c]onformance\.sh'; do sleep 30; done

log "simulated reads"; "$ROOT/bench/sim_wg_eval.sh" "$R/sim_wg" > "$R/sim_wg.log" 2>&1
for set in sup20 hac20 r9_20; do
  case $set in
    sup20) reads="$ROOT/data/reads/hg002_PAW70337_chr20.fq.gz"; model=r1041_e82_400bps_sup_v500 ;;
    hac20) reads="$G/hg002_PAW70337_hac_chr20.fq.gz"; model=r1041_e82_400bps_hac_v500 ;;
    r9_20) reads="$G/hg002_r941_sup_chr20.fq.gz"; model=r941_prom_sup_g5014 ;;
  esac
  for a in minimap2 bonefish; do
    label="v012_${set}_$a"; al=""; [[ $a == minimap2 ]] && al=minimap2
    [[ -s "$R/accuracy/$label/sv_accuracy.tsv" ]] && continue
    log "accuracy $label"
    keep=""; [[ $set == sup20 ]] && keep=1
    KEEP_BAM=$keep ALIGNER=$al READS="$reads" CLAIR3_MODEL=$model "$ROOT/bench/eval_accuracy.sh" "$label" "" > "$R/accuracy_$label.log" 2>&1
  done
done
log "phasing"; "$ROOT/bench/phasing.sh" "$R/phasing" v012_sup20_minimap2 v012_sup20_bonefish > "$R/phasing.log" 2>&1
log "CMRG and tandem repeats"; "$ROOT/bench/cmrg_tr_eval.sh" "$R/cmrg_tr" > "$R/cmrg_tr.log" 2>&1
log "determinism"; "$ROOT/bench/determinism.sh" "$R/determinism" > "$R/determinism.log" 2>&1
log "large reference"; "$ROOT/bench/large_ref.sh" "$R/large_ref" > "$R/large_ref.log" 2>&1
log "scaling (waiting for an idle machine)"
while [[ $(awk '{print int($1)}' /proc/loadavg) -gt 1 ]]; do sleep 30; done
"$ROOT/bench/scaling.sh" "$R/scaling" > "$R/scaling.log" 2>&1
log "done"
