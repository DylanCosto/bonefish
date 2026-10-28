#!/usr/bin/env bash
# Whole-genome HG002 HiFi (PacBio Revio 2022Q4 rep1, ~30x): bonefish (default dna preset) vs minimap2
# 2.31 map-hifi (lr:hq optional: ALIGNERS=...mm2_lrhq) on GRCh38 no-alt. Per aligner: aligner process time and whole-pipeline
# wall (aligner | samtools sort), Clair3 hifi_revio vs GIAB v4.2.1 (SNPs, indels), Sniffles2 + Truvari
# vs GIAB v5.0q SVs, primary placement vs minimap2 map-hifi. bonefish_hifi: -x hifi (ALIGNERS=bonefish_hifi).
set -uo pipefail
R=${BENCH}; H=$R/hifi-data; W=$R/hifi-work; O=$H/wgs_eval; mkdir -p $O
VE=$HOME/micromamba/envs/vareval/bin; export PATH=$VE:$HOME/micromamba/envs/ontaln/bin:$PATH
BF=${BF:-$R/bin/bonefish-v030}; BFH=${BFH:-$R/bin/bonefish-hifi-eval3}; T=16; MODEL=$VE/models/hifi_revio
REF=$W/ref/GRCh38_no_alt.fa; SDF=$W/ref/GRCh38.sdf; Q=$H/wgs/HG002.hifi.fq.gz
BIDX=$H/index/GRCh38.bonefish.raw; MIDX=$H/index/GRCh38.mm2.k19w19.mmi
TV=$H/truth/wgs/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz; TB=$H/truth/wgs/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed
SVB=$H/truth/wgs/HG002_GRCh38_v5.0q_stvar.benchmark.bed; SVT=$O/truth_sv.vcf.gz
ALIGNERS=${ALIGNERS:-"bonefish mm2_hifi"}
log() { echo "[$(date '+%m-%d %H:%M:%S')] $*" | tee -a $O/eval.log; }
for k in snps indels; do
  [[ -s $O/truth_$k.vcf.gz ]] || { bcftools view -v $k -Oz -o $O/truth_$k.vcf.gz $TV && tabix -p vcf $O/truth_$k.vcf.gz; }
done
[[ -s $SVT ]] || { bcftools norm -m -any -Oz -o $SVT $H/truth/wgs/HG002_GRCh38_v5.0q_stvar.vcf.gz && tabix -p vcf $SVT; }
for a in $ALIGNERS; do
  D=$O/$a; mkdir -p $D; BAM=$D/aln.bam
  if [[ ! -s $BAM.bai ]]; then
    case $a in
      mm2_hifi) IDX=$MIDX; cmd="minimap2 -ax map-hifi -y -t $T $MIDX $Q" ;;
      mm2_lrhq) IDX=$MIDX; cmd="minimap2 -ax lr:hq -y -t $T $MIDX $Q" ;;
      bonefish) IDX=$BIDX; cmd="$BF -a -y -t $T $BIDX $Q" ;;
      bonefish_hifi) IDX=$BIDX; cmd="$BFH -x hifi -a -y -t $T $BIDX $Q" ;;
    esac
    cat $IDX > /dev/null
    log "$a: aligning"
    /usr/bin/time -f "%e %U %S %M" -o $D/pipeline.time bash -c "set -o pipefail; /usr/bin/time -f '%e %U %S %M' -o $D/align.time $cmd 2> $D/align.err | $H/bin/hifiwgs-sort sort -@ 4 -m 1G -T $D/sorttmp -o $BAM - 2> $D/sort.err"
    samtools index -@ 8 $BAM
    read -r w u sy m < $D/align.time; read -r pw _ < $D/pipeline.time
    log "$a: align wall=$w cpu=$(echo "$u+$sy"|bc) rss=$((m/1024))MB; reads->sorted BAM wall=$pw"
  fi
  if [[ ! -s $D/clair3/merge_output.vcf.gz ]]; then
    rm -rf $D/clair3; log "$a: clair3"
    /usr/bin/time -f "%e %U %S %M" -o $D/clair3.time run_clair3.sh --bam_fn=$BAM --ref_fn=$REF --threads=$T \
      --platform=hifi --model_path=$MODEL --output=$D/clair3 \
      --use_longphase_for_intermediate_phasing=True --longphase=$VE/longphase > $D/clair3.log 2>&1
    log "$a: clair3 $(cut -d' ' -f1 $D/clair3.time) s"
    rm -rf $D/clair3/tmp
  fi
  if [[ ! -s $D/accuracy.tsv ]]; then
    echo -e "type\ttp\tfp\tfn\tprecision\trecall\tf1" > $D/accuracy.tsv
    for k in snps indels; do
      bcftools view -v $k -f PASS -Oz -o $D/calls_$k.vcf.gz $D/clair3/merge_output.vcf.gz; tabix -f -p vcf $D/calls_$k.vcf.gz
      rm -rf $D/vcfeval_$k
      rtg vcfeval -b $O/truth_$k.vcf.gz -c $D/calls_$k.vcf.gz -e $TB -t $SDF -o $D/vcfeval_$k --threads $T > /dev/null
      awk -v k=$k '$1=="None"{printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n", k, $2, $4, $5, $6, $7, $8}' $D/vcfeval_$k/summary.txt >> $D/accuracy.tsv
    done
    log "$a accuracy: $(tail -n +2 $D/accuracy.tsv | awk '{printf "%s F1 %s (FP %s FN %s); ", $1, $7, $3, $4}')"
  fi
  if [[ ! -s $D/sv_accuracy.tsv ]]; then
    sniffles --input $BAM --vcf $D/sv.vcf --reference $REF --threads $T > $D/sniffles.log 2>&1
    bcftools sort -Oz -o $D/sv.vcf.gz $D/sv.vcf 2>/dev/null && tabix -f -p vcf $D/sv.vcf.gz && rm -f $D/sv.vcf
    echo -e "matching\ttp\tfp\tfn\tprecision\trecall\tf1" > $D/sv_accuracy.tsv
    for mode in strict multi; do
      extra=(); [[ $mode == multi ]] && extra=(--pick multi)
      rm -rf $D/truvari_$mode
      truvari bench -b $SVT -c $D/sv.vcf.gz --includebed $SVB --passonly --pctseq 0 "${extra[@]}" -o $D/truvari_$mode > $D/truvari_$mode.log 2>&1
      python3 -c "import json; s=json.load(open('$D/truvari_$mode/summary.json')); print('\t'.join(['$mode',str(s['TP-base']),str(s['FP']),str(s['FN'])]+[f\"{s[k]:.4f}\" for k in ('precision','recall','f1')]))" >> $D/sv_accuracy.tsv
    done
    log "$a SV: $(tail -n +2 $D/sv_accuracy.tsv | awk '{printf "%s F1 %s (TP %s FP %s FN %s); ", $1, $7, $2, $3, $4}')"
  fi
done
for a in bonefish bonefish_hifi mm2_lrhq; do
  if [[ -s $O/mm2_hifi/aln.bam && -s $O/$a/aln.bam && ! -s $O/$a/compare_vs_mm2hifi.txt ]]; then
    python3 $R/bench/compare_primary.py $O/mm2_hifi/aln.bam $O/$a/aln.bam > $O/$a/compare_vs_mm2hifi.txt 2>/dev/null
    log "$a vs mm2_hifi: $(grep -E 'identical|same locus|moved|unmapped' $O/$a/compare_vs_mm2hifi.txt | tr -s ' ' | tr '\n' ' ')"
  fi
done
log DONE
