#!/bin/bash

set -e

function usage {
  echo "WORKDIR= BAM_DIR= run_hapcaller.sh"
  exit 1
}

function print_abs_path {
  echo $(cd $(dirname $1); pwd)/$(basename $1)
}

[ -z "$WORKDIR" -o -z "$BAM_DIR" ] && usage

mkdir -pv $WORKDIR

WORKDIR=$(print_abs_path $WORKDIR)
BAM_DIR=$(print_abs_path $BAM_DIR)

reference=~/lake/users/marghoob/GATK-bundle-hg19/ucsc.hg19.fa
dbsnp=~/lake/users/marghoob/GATK-bundle-hg19/dbsnp_137.hg19.vcf
GATK_JAR=~/lake/opt/gatk-2.7-2-g6bda569/GenomeAnalysisTK.jar

if [ -z "$CHR_LIST" ]; then
  CHR_LIST=`awk '{print $1}' $reference.fai`
fi

mkdir -pv $WORKDIR/beds
mkdir -pv $WORKDIR/bams
mkdir -pv $WORKDIR/vcfs
mkdir -pv $WORKDIR/logs

START=$(date +%s)

echo "Running haplotype caller on $BAM_DIR"

for chr in $CHR_LIST; do
  [ ! -s "$BAM_DIR/$chr.bam" ] && continue
  (cd $WORKDIR/bams && ln -sf $BAM_DIR/$chr.bam && ln -sf $BAM_DIR/$chr.bai) 
  awk -v chr=$chr '{if ($1 == chr) {print chr "\t0\t" $2; exit}}' $reference.fai > $WORKDIR/beds/$chr.bed
  echo "java -Xmx5g -Xms5g -jar $GATK_JAR -T HaplotypeCaller -L $WORKDIR/beds/$chr.bed --input_file $WORKDIR/bams/$chr.bam --out $WORKDIR/vcfs/$chr.vcf -R $reference --dbsnp $dbsnp -stand_emit_conf 0.1 -A AlleleBalance -A Coverage -A MappingQualityZero &>$WORKDIR/logs/$chr.log"
done | xargs -I CMD --max-procs=12 bash -c CMD

vcf_list=
for chr in $CHR_LIST; do
  [ ! -s "$BAM_DIR/$chr.bam" ] && continue
  vcf_list="$vcf_list $WORKDIR/vcfs/$chr.vcf"
done

echo "Will merge $vcf_list"
# now merge the vcfs
cat $vcf_list | awk 'BEGIN {header_seen = 0} /^#/ {if (header_seen == 0) print $0} !/^#/ {header_seen = 1; print $0}' > $WORKDIR/vcfs/merged.vcf

# Split into SNPs and INDELS

echo "Will split into SNPs and INDELs"
for vartype in SNP INDEL; do
  java -jar $GATK_JAR -T SelectVariants -selectType $vartype -V $WORKDIR/vcfs/merged.vcf -o $WORKDIR/vcfs/merged.$vartype.vcf -R $reference &>$WORKDIR/logs/SelectVariants.$vartype.log &
done
wait

# VQSR
echo "Will run VariantRecalibrator now"
java -jar $GATK_JAR -T VariantRecalibrator -input $WORKDIR/vcfs/merged.SNP.vcf -resource:hapmap,VCF,known=false,training=true,truth=true,prior=15.0 ~/lake/users/marghoob/GATK-bundle-hg19/hapmap_3.3.hg19.vcf -resource:omni,VCF,known=false,training=true,truth=true,prior=12.0 ~/lake/users/marghoob/GATK-bundle-hg19/1000G_omni2.5.hg19.vcf -resource:1000G,VCF,known=false,training=true,truth=false,prior=10.0 ~/lake/users/marghoob/GATK-bundle-hg19/1000G_phase1.snps.high_confidence.hg19.vcf \
 -resource:dbsnp,VCF,known=true,training=false,truth=false,prior=2.0 ~/lake/users/marghoob/GATK-bundle-hg19/dbsnp_137.hg19.vcf -an QD -an MQRankSum -an ReadPosRankSum -an FS -an DP -mode SNP --tranches_file $WORKDIR/vcfs/SNP.tranches --recal_file $WORKDIR/vcfs/SNP.recal.vcf -R $reference &>$WORKDIR/logs/VariantRecalibrator.SNP.log &
java -jar $GATK_JAR -T VariantRecalibrator -input $WORKDIR/vcfs/merged.INDEL.vcf -resource:mills,VCF,known=false,training=true,truth=true,prior=12.0 ~/lake/users/marghoob/GATK-bundle-hg19/Mills_and_1000G_gold_standard.indels.hg19.vcf -resource:dbsnp,VCF,known=true,training=false,truth=false,prior=2.0 ~/lake/users/marghoob/GATK-bundle-hg19/dbsnp_137.hg19.vcf -an DP -an FS -an MQRankSum -an ReadPosRankSum -mode INDEL --maxGaussians 4 --tranches_file $WORKDIR/vcfs/INDEL.tranches --recal_file $WORKDIR/vcfs/INDEL.recal.vcf -R $reference &>$WORKDIR/logs/VariantRecalibrator.INDEL.log &
wait

echo "Will run ApplyRecalibration now"
for vartype in SNP INDEL
do
 java -jar $GATK_JAR -T ApplyRecalibration -mode $vartype -input $WORKDIR/vcfs/merged.$vartype.vcf --ts_filter_level 99.9 --tranches_file $WORKDIR/vcfs/$vartype.tranches --recal_file $WORKDIR/vcfs/$vartype.recal.vcf -R $reference -o $WORKDIR/vcfs/$vartype.recalibrated.vcf &>$WORKDIR/logs/ApplyRecalibration.$vartype.log &
done
wait

echo "Done"
END=$(date +%s)
DIFF=$(( $END - $START ))
echo "It took $DIFF seconds" > $WORKDIR/logs/time.log

