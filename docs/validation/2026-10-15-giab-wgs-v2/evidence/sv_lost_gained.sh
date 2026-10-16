#!/usr/bin/env bash
# Truth SVs found (Truvari strict TP) by one run and missed by another, by type, size and tandem-repeat overlap.
# Usage: bench/sv_lost_gained.sh OLD_TRUVARI_DIR NEW_TRUVARI_DIR TR_TIER1_BED
set -euo pipefail
export PATH=$HOME/micromamba/envs/vareval/bin:$PATH
A=$1; B=$2; TR=$3; T=$(mktemp -d); trap 'rm -rf $T' EXIT
bcftools isec -c none -n=2 -w1 -p $T/lost "$A/tp-base.vcf.gz" "$B/fn.vcf.gz" > /dev/null 2>&1
bcftools isec -c none -n=2 -w1 -p $T/gained "$A/fn.vcf.gz" "$B/tp-base.vcf.gz" > /dev/null 2>&1
for k in lost gained; do
  f=$T/$k/0000.vcf
  echo -e "## $k (truth SVs, new run vs old)\ttotal $(grep -vc '^#' $f)\tin Tier 1 tandem repeats $(bcftools view -H -T "$TR" $f | wc -l)"
  grep -v '^#' $f | awk '{l=length($5)-length($4); if(l<0)l=-l; t=(length($5)>length($4))?"INS":"DEL";
    b=(l<100)?"50-99":(l<300)?"100-299":(l<1000)?"300-999":">=1000"; print t"\t"b}' | sort | uniq -c | sort -rn
done
