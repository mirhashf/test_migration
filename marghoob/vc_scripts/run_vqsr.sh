#!/bin/bash

set -ex

REFERENCE=~/lake/users/marghoob/GATK-bundle-hg19/ucsc.hg19.fa
GATK_JAR=~/lake/opt/gatk-2.7-2-g6bda569/GenomeAnalysisTK.jar

java -jar $GATK_JAR -T VariantRecalibrator -input merged.SNP.vcf.gz -resource:hapmap,VCF,known=false,training=true,truth=true,prior=15.0 ~/lake/users/marghoob/GATK-bundle-hg19/hapmap_3.3.hg19.vcf -resource:omni,VCF,known=false,training=true,truth=true,prior=12.0 ~/lake/users/marghoob/GATK-bundle-hg19/1000G_omni2.5.hg19.vcf -resource:1000G,VCF,known=false,training=true,truth=false,prior=10.0 ~/lake/users/marghoob/GATK-bundle-hg19/1000G_phase1.snps.high_confidence.hg19.vcf -resource:dbsnp,VCF,known=true,training=false,truth=false,prior=2.0 ~/lake/users/marghoob/GATK-bundle-hg19/dbsnp_137.hg19.vcf -an QD -an MQRankSum -an ReadPosRankSum -an FS -an DP -mode SNP --tranches_file SNP.tranches --recal_file SNP.recal.vcf -R $REFERENCE &>VariantRecalibrator.SNP.log &
java -jar $GATK_JAR -T VariantRecalibrator -input merged.INDEL.vcf.gz -resource:mills,VCF,known=false,training=true,truth=true,prior=12.0 ~/lake/users/marghoob/GATK-bundle-hg19/Mills_and_1000G_gold_standard.indels.hg19.vcf -resource:dbsnp,VCF,known=true,training=false,truth=false,prior=2.0 ~/lake/users/marghoob/GATK-bundle-hg19/dbsnp_137.hg19.vcf -an DP -an FS -an MQRankSum -an ReadPosRankSum -mode INDEL --maxGaussians 4 --tranches_file INDEL.tranches --recal_file INDEL.recal.vcf -R $REFERENCE &>VariantRecalibrator.INDEL.log &
wait

for vartype in SNP INDEL
do
java -jar $GATK_JAR -T ApplyRecalibration -mode $vartype -input merged.$vartype.vcf.gz --ts_filter_level 99.9 --tranches_file $vartype.tranches --recal_file $vartype.recal.vcf -R $REFERENCE -o $vartype.recalibrated.vcf &>ApplyRecalibration.$vartype.log &
done
wait

bgzip -f SNP.recalibrated.vcf && tabix -f SNP.recalibrated.vcf.gz
bgzip -f INDEL.recalibrated.vcf && tabix -f INDEL.recalibrated.vcf.gz
