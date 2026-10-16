#!/usr/bin/env bash
# Are changed indel calls between two chr20 runs tied to changed alignments? For new FPs and lost calls
# of run B vs run A, and for random truth indels, the share with any primary-alignment change
# (name, flag, pos, MAPQ, CIGAR) within 300 bp. Usage: bench/indel_churn.sh A_DIR B_DIR A_BAM B_BAM
set -euo pipefail
export PATH=$HOME/micromamba/envs/vareval/bin:$PATH
A=$1; B=$2; ABAM=$3; BBAM=$4; T=$(mktemp -d); trap 'rm -rf $T' EXIT
bcftools isec -c none -p $T/fp "$B/vcfeval_indels/fp.vcf.gz" "$A/calls_indels.vcf.gz" > /dev/null 2>&1
bcftools isec -c none -p $T/fn "$B/vcfeval_indels/fn.vcf.gz" "$A/calls_indels.vcf.gz" > /dev/null 2>&1
bcftools view -H "$B/vcfeval_indels/tp-baseline.vcf.gz" | shuf -n 300 --random-source=<(yes) > $T/rand.vcf
changed() { local n=0 c=0; for p in $(grep -v '^#' "$1" | awk '{print $1":"($2>300?$2-300:1)"-"$2+300}'); do
  n=$((n+1)); [[ "$(samtools view -F 0x900 "$ABAM" $p | cut -f1-6 | sort | md5sum)" != \
                 "$(samtools view -F 0x900 "$BBAM" $p | cut -f1-6 | sort | md5sum)" ]] && c=$((c+1)); done
  echo -e "$2\t$n\t$c"; }
echo -e "set\tsites\twith_alignment_change_300bp"
changed $T/fp/0000.vcf "new FPs (B false positive, not called by A)"
changed $T/fn/0002.vcf "lost calls (B false negative, called by A)"
changed $T/rand.vcf "random truth indels"
