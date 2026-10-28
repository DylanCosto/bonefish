#!/usr/bin/env bash
# Long-gap cluster realignment (-x hifi) on ONT reads: HG002 R10.4.1 sup v5.0.0 (PAW70337), reads
# with a primary alignment on chr20, against the whole GRCh38 genome. bonefish dna (v0.3.0) vs the
# realignment at identity thresholds 0.9 (the hifi preset), 0.8 and 0 (gaps joined, never rewritten
# as deletion + insertion), and minimap2 2.31 lr:hq. Per run: aligner time, Sniffles2 + Truvari vs
# GIAB v5.0q stvar (chr20), Clair3 r1041_e82_400bps_sup_v500 vs GIAB v4.2.1 (chr20).
set -uo pipefail
R=${BENCH}; H=$R/hifi-data; W=$R/ont-data/realign_eval; mkdir -p $W
VE=$HOME/micromamba/envs/vareval/bin; export PATH=$VE:$HOME/micromamba/envs/ontaln/bin:$PATH
BF=$R/bin/bonefish-v030; BFX=$R/bin/bonefish-hifi-eval3; T=16; MODEL=$VE/models/r1041_e82_400bps_sup_v500
REF=$R/hifi-work/ref/GRCh38_no_alt.fa; SDF=$R/hifi-work/ref/GRCh38.sdf; Q=$R/ont-data/reads/hg002_PAW70337_chr20.fq.gz
BIDX=$H/index/GRCh38.bonefish.raw; MIDX=$H/index/GRCh38.mm2.k19w19.mmi
ALIGNERS=${ALIGNERS:-"bonefish_dna bonefish_rg90 bonefish_rg80 bonefish_rg0 mm2_lrhq"}
log() { echo "[$(date '+%m-%d %H:%M:%S')] $*" | tee -a $W/eval.log; }
SVT=$W/truth_sv.vcf.gz
[[ -s $SVT ]] || { bcftools norm -m -any -Oz -o $SVT $H/truth/HG002_chr20_v5.0q_stvar.vcf.gz && tabix -p vcf $SVT; }
for k in snps indels; do
  [[ -s $W/truth_$k.vcf.gz ]] || { bcftools view -v $k -Oz -o $W/truth_$k.vcf.gz $H/truth/HG002_chr20_v4.2.1.vcf.gz && tabix -p vcf $W/truth_$k.vcf.gz; }
done
for a in $ALIGNERS; do
  D=$W/$a; mkdir -p $D; BAM=$D/aln.bam
  if [[ ! -s $BAM.bai ]]; then
    env=()
    case $a in
      bonefish_dna) IDX=$BIDX; cmd=($BF -a -y -t $T $BIDX $Q) ;;
      bonefish_rg90) IDX=$BIDX; cmd=($BFX -x hifi -a -y -t $T $BIDX $Q) ;;
      bonefish_rg80) IDX=$BIDX; env=(BONEFISH_EXP_REALIGN_MIN_ID=0.8); cmd=($BFX -x hifi -a -y -t $T $BIDX $Q) ;;
      bonefish_rg0) IDX=$BIDX; env=(BONEFISH_EXP_REALIGN_MIN_ID=0); cmd=($BFX -x hifi -a -y -t $T $BIDX $Q) ;;
      mm2_lrhq) IDX=$MIDX; cmd=(minimap2 -ax lr:hq -y -t $T $MIDX $Q) ;;
    esac
    cat $IDX > /dev/null
    env "${env[@]}" /usr/bin/time -f "%e %U %S %M" -o $D/align.time "${cmd[@]}" 2> $D/align.err \
      | samtools sort -@ 4 -m 1G -T $D/sorttmp -o $BAM - 2> $D/sort.err
    samtools index $BAM
    read -r w u sy m < $D/align.time; log "$a align wall=$w cpu=$(echo "$u+$sy"|bc) rss=$((m/1024))MB"
  fi
  if [[ ! -s $D/sv_accuracy.tsv ]]; then
    sniffles --input $BAM --vcf $D/sv.vcf --reference $REF --threads $T > $D/sniffles.log 2>&1
    bcftools sort -Oz -o $D/sv.vcf.gz $D/sv.vcf 2>/dev/null && tabix -f -p vcf $D/sv.vcf.gz && rm -f $D/sv.vcf
    echo -e "matching\ttp\tfp\tfn\tprecision\trecall\tf1" > $D/sv_accuracy.tsv
    for mode in strict multi; do
      extra=(); [[ $mode == multi ]] && extra=(--pick multi)
      rm -rf $D/truvari_$mode
      truvari bench -b $SVT -c $D/sv.vcf.gz --includebed $H/truth/HG002_chr20_v5.0q_stvar.bed --passonly --pctseq 0 "${extra[@]}" -o $D/truvari_$mode > $D/truvari_$mode.log 2>&1
      python3 -c "import json,sys; d=json.load(open(sys.argv[1])); print('\t'.join(str(x) for x in [sys.argv[2], d['TP-comp'], d['FP'], d['FN'], round(d['precision'],4), round(d['recall'],4), round(d['f1'],4)]))" $D/truvari_$mode/summary.json $mode >> $D/sv_accuracy.tsv
    done
    log "$a SV: $(tail -n +2 $D/sv_accuracy.tsv | awk '{printf "%s F1 %s (TP %s FP %s FN %s); ", $1, $7, $2, $3, $4}')"
  fi
  if [[ $a != bonefish_rg80 && $a != bonefish_rg0 && ! -s $D/accuracy.tsv ]]; then
    rm -rf $D/clair3
    /usr/bin/time -f "%e %U %S %M" -o $D/clair3.time run_clair3.sh --bam_fn=$BAM --ref_fn=$REF --threads=$T \
      --platform=ont --model_path=$MODEL --output=$D/clair3 --ctg_name=chr20 \
      --use_longphase_for_intermediate_phasing=True --longphase=$VE/longphase > $D/clair3.log 2>&1
    echo -e "type\ttp\tfp\tfn\tprecision\trecall\tf1" > $D/accuracy.tsv
    for k in snps indels; do
      bcftools view -v $k -f PASS -Oz -o $D/calls_$k.vcf.gz $D/clair3/merge_output.vcf.gz; tabix -f -p vcf $D/calls_$k.vcf.gz
      rm -rf $D/vcfeval_$k
      rtg vcfeval -b $W/truth_$k.vcf.gz -c $D/calls_$k.vcf.gz -e $H/truth/HG002_chr20_v4.2.1.bed -t $SDF -o $D/vcfeval_$k --threads $T > /dev/null
      awk -v k=$k '$1=="None"{printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n", k, $2, $4, $5, $6, $7, $8}' $D/vcfeval_$k/summary.txt >> $D/accuracy.tsv
    done
    log "$a accuracy: $(tail -n +2 $D/accuracy.tsv | awk '{printf "%s F1 %s (FP %s FN %s); ", $1, $7, $3, $4}')"
  fi
done
log DONE
