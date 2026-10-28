#!/usr/bin/env bash
# chr20 HiFi (PacBio Revio, HG002/HG003/HG004) against the whole GRCh38 genome: minimap2 2.31
# map-hifi and lr:hq vs bonefish (default dna preset; bonefish_hifi: -x hifi). Per run: aligner time (SAM piped to samtools
# sort, aligner process only), Clair3 hifi_revio small variants vs GIAB v4.2.1 (chr20), and for
# HG002 Sniffles2 + Truvari vs GIAB v5.0q stvar (chr20); primary placement vs minimap2 map-hifi.
set -uo pipefail
R=${BENCH}; H=$R/hifi-data; W=$R/hifi-work; mkdir -p $W/ref
VE=$HOME/micromamba/envs/vareval/bin; export PATH=$VE:$HOME/micromamba/envs/ontaln/bin:$PATH
BF=${BF:-$R/bin/bonefish-v030}; BFH=${BFH:-$R/bin/bonefish-hifi-eval3}; T=16; MODEL=$VE/models/hifi_revio
SAMPLES=${SAMPLES:-"HG002 HG003 HG004"}; ALIGNERS=${ALIGNERS:-"mm2_hifi mm2_lrhq bonefish"}
REF=$W/ref/GRCh38_no_alt.fa; SDF=$W/ref/GRCh38.sdf
log() { echo "[$(date +%H:%M:%S)] $*" | tee -a $W/eval.log; }
[[ -s $REF.fai ]] || { cp $H/ref/GRCh38_no_alt.fa $H/ref/GRCh38_no_alt.fa.fai $W/ref/; }
[[ -d $SDF ]] || rtg format -o $SDF $REF > /dev/null
BIDX=$H/index/GRCh38.bonefish.raw; MIDX=$H/index/GRCh38.mm2.k19w19.mmi
[[ -s $BIDX ]] || { /usr/bin/time -f "%e s %M KB" -o $H/index/bonefish.idx.time $BF -t $T -d $BIDX --idx-raw $REF 2> $H/index/bonefish.idx.log; }
[[ -s $MIDX ]] || { /usr/bin/time -f "%e s %M KB" -o $H/index/mm2.idx.time minimap2 -x map-hifi -t $T -d $MIDX $REF 2> $H/index/mm2.idx.log; }
SVT=$W/truth_HG002_sv.vcf.gz
[[ -s $SVT ]] || { bcftools norm -m -any -Oz -o $SVT $H/truth/HG002_chr20_v5.0q_stvar.vcf.gz && tabix -p vcf $SVT; }
for s in $SAMPLES; do
  for k in snps indels; do
    [[ -s $W/truth_${s}_$k.vcf.gz ]] || { bcftools view -v $k -Oz -o $W/truth_${s}_$k.vcf.gz $H/truth/${s}_chr20_v4.2.1.vcf.gz && tabix -p vcf $W/truth_${s}_$k.vcf.gz; }
  done
  for a in $ALIGNERS; do
    D=$W/$s/$a; mkdir -p $D; BAM=$D/aln.bam; Q=$H/reads/$s.chr20.fq.gz
    if [[ ! -s $BAM.bai ]]; then
      case $a in
        mm2_hifi) cat $MIDX > /dev/null; cmd=(minimap2 -ax map-hifi -y -t $T $MIDX $Q) ;;
        mm2_lrhq) cat $MIDX > /dev/null; cmd=(minimap2 -ax lr:hq -y -t $T $MIDX $Q) ;;
        bonefish) cat $BIDX > /dev/null; cmd=($BF -a -y -t $T $BIDX $Q) ;;
        bonefish_hifi) cat $BIDX > /dev/null; cmd=($BFH -x hifi -a -y -t $T $BIDX $Q) ;;
      esac
      /usr/bin/time -f "%e %U %S %M" -o $D/align.time "${cmd[@]}" 2> $D/align.err | samtools sort -@ 4 -m 1500M -o $BAM - 2>/dev/null
      samtools index $BAM
      read -r w u sy m < $D/align.time; log "$s $a align wall=$w cpu=$(echo "$u+$sy"|bc) rss=$((m/1024))MB"
    fi
    if [[ ! -s $D/clair3/merge_output.vcf.gz ]]; then
      rm -rf $D/clair3
      /usr/bin/time -f "%e %U %S %M" -o $D/clair3.time run_clair3.sh --bam_fn=$BAM --ref_fn=$REF --threads=$T \
        --platform=hifi --model_path=$MODEL --output=$D/clair3 --ctg_name=chr20 \
        --use_longphase_for_intermediate_phasing=True --longphase=$VE/longphase > $D/clair3.log 2>&1
      log "$s $a clair3 $(cut -d' ' -f1 $D/clair3.time) s"
    fi
    if [[ ! -s $D/accuracy.tsv ]]; then
      echo -e "type\ttp\tfp\tfn\tprecision\trecall\tf1" > $D/accuracy.tsv
      for k in snps indels; do
        bcftools view -v $k -f PASS -Oz -o $D/calls_$k.vcf.gz $D/clair3/merge_output.vcf.gz; tabix -f -p vcf $D/calls_$k.vcf.gz
        rm -rf $D/vcfeval_$k
        rtg vcfeval -b $W/truth_${s}_$k.vcf.gz -c $D/calls_$k.vcf.gz -e $H/truth/${s}_chr20_v4.2.1.bed -t $SDF -o $D/vcfeval_$k --threads $T > /dev/null
        awk -v k=$k '$1=="None"{printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n", k, $2, $4, $5, $6, $7, $8}' $D/vcfeval_$k/summary.txt >> $D/accuracy.tsv
      done
      log "$s $a accuracy: $(tail -n +2 $D/accuracy.tsv | awk '{printf "%s F1 %s (FP %s FN %s); ", $1, $7, $3, $4}')"
    fi
    if [[ $s == HG002 && ! -s $D/sv_accuracy.tsv ]]; then
      sniffles --input $BAM --vcf $D/sv.vcf --reference $REF --threads $T > $D/sniffles.log 2>&1
      bcftools sort -Oz -o $D/sv.vcf.gz $D/sv.vcf 2>/dev/null && tabix -f -p vcf $D/sv.vcf.gz && rm -f $D/sv.vcf
      echo -e "matching\ttp\tfp\tfn\tprecision\trecall\tf1" > $D/sv_accuracy.tsv
      for mode in strict multi; do
        extra=(); [[ $mode == multi ]] && extra=(--pick multi)
        rm -rf $D/truvari_$mode
        truvari bench -b $SVT -c $D/sv.vcf.gz --includebed $H/truth/HG002_chr20_v5.0q_stvar.bed --passonly --pctseq 0 "${extra[@]}" -o $D/truvari_$mode > $D/truvari_$mode.log 2>&1
        python3 -c "import json; s=json.load(open('$D/truvari_$mode/summary.json')); print('\t'.join(['$mode',str(s['TP-base']),str(s['FP']),str(s['FN'])]+[f\"{s[k]:.4f}\" for k in ('precision','recall','f1')]))" >> $D/sv_accuracy.tsv
      done
      log "$s $a SV: $(tail -n +2 $D/sv_accuracy.tsv | awk '{printf "%s F1 %s (TP %s FP %s FN %s); ", $1, $7, $2, $3, $4}')"
    fi
    if [[ $a != mm2_hifi && -s $W/$s/mm2_hifi/aln.bam && ! -s $D/compare_vs_mm2hifi.txt ]]; then
      python3 $R/bench/compare_primary.py $W/$s/mm2_hifi/aln.bam $BAM > $D/compare_vs_mm2hifi.txt 2>/dev/null
      log "$s $a vs mm2_hifi: $(grep -E 'identical|unmapped' $D/compare_vs_mm2hifi.txt | tr -s ' ' | tr '\n' ' ')"
    fi
  done
done
log DONE
