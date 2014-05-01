#!/bin/bash

set -ex

myname=`basename $0`
function usage {
  echo "CHR_LIST= $myname <outdir>"
  echo "CHR_LIST is optional. If it is empty, all chromosomes are generated."
  exit 1
}

function print_abs_path {
  echo $(cd $(dirname $1); pwd)/$(basename $1)
}

LAKE=/net/kodiak/volumes/lake/shared

NIST=$LAKE/users/marghoob/NIST/NISThighConf
GATKVCF=$LAKE/users/marghoob/GATK-bundle-hg19/CEUTrio.HiSeq.WGS.b37.bestPractices.phased.hg19.vcf.gz
DELETIONS_VCF=$LAKE/users/marghoob/homes.gersteinlab.org/people/aabyzov/outbox/private/NA12878_variants/NA12878.2010_06.and.fosmid.deletions.phased.vcf

OUTDIR=$1
[ -z "$OUTDIR" ] && usage

mkdir -pv $OUTDIR

OUTDIR=$(print_abs_path $OUTDIR)

VCFDIR=$OUTDIR/vcfs
GENOMEDIR=$OUTDIR/genome
LOGDIR=$OUTDIR/log
mkdir -pv $OUTDIR $VCFDIR $GENOMEDIR $LOGDIR

# Tool vars
DIR="$( cd "$( dirname "$0" )" && pwd )"
LIFTVCF=$DIR/../lift_vcf.sh
GATK_JAR=$LAKE/opt/CancerAnalysisPackage-2013.2-18-g8207e53/GenomeAnalysisTK.jar
VCF2DIPLOID=$LAKE/users/marghoob/vcf2diploid/vcf2diploid.jar
IGVTOOLS=$LAKE/users/marghoob/IGVTools/igvtools.jar
DWGSIM=$LAKE/opt/dwgsim/dwgsim
BEDTOOLS_DIR=$LAKE/opt/bedtools-2.17.0/bin/
SAMTOOLS=$LAKE/opt/samtools/samtools
export PATH=$LAKE/opt/tabix-0.2.6:$LAKE/opt/vcftools_0.1.11/bin:$PATH
export PERL5LIB=$LAKE/opt/vcftools_0.1.11/perl

REFERENCE=$LAKE/users/marghoob/GATK-bundle-hg19/ucsc.hg19.fa

mkdir -pv $VCFDIR/tmp
cd $VCFDIR/tmp

# Prepare the right VCF for deletions
echo "Lifting the deletion SVs from b36 to b37 after removing records with incorrect filter and removing the END field in INFO column"
awk '/^#/ { str = $0; sub(/-1/, "1", str); print str} !/^#/ { str = $0; sub(/END=[0-9]+;/, "", str); if ($7 != "0") print str}' $DELETIONS_VCF > deletions.b36.vcf
$LIFTVCF b36 b37 deletions.b36.vcf deletions.b37.vcf &>$LOGDIR/deletions.b36.to.b37.log
echo "Lifting the deletion SVs from b37 to hg19"
$LIFTVCF b37 hg19 deletions.b37.vcf deletions.hg19.vcf &>$LOGDIR/deletions.b37.to.hg19.log
bgzip -f deletions.hg19.vcf && tabix -f deletions.hg19.vcf.gz

cd ..

# Extract only the relevant chromosomes from the VCFs
if [ -n "$CHR_LIST" ]; then
  (tabix -h tmp/deletions.hg19.vcf.gz $CHR_LIST | bgzip > deletions.vcf.gz && tabix -f deletions.vcf.gz)
  (tabix -h $GATKVCF $CHR_LIST | bgzip > GATK.vcf.gz && tabix -f GATK.vcf.gz)
  for var in SNP INDEL; do
    (tabix -h $NIST.$var.hg19.annotated.vcf.gz $CHR_LIST | bgzip > NIST.$var.vcf.gz && tabix -f NIST.$var.vcf.gz)
  done
else
  ln -sf tmp/deletions.hg19.vcf.gz deletions.vcf.gz && ln -sf tmp/deletions.hg19.vcf.gz.tbi deletions.vcf.gz.tbi
  ln -sf $GATKVCF GATK.vcf.gz && ln -sf $GATKVCF.tbi GATK.vcf.gz.tbi
  for var in SNP INDEL; do
    ln -sf $NIST.$var.hg19.annotated.vcf.gz NIST.$var.vcf.gz && ln -sf $NIST.$var.hg19.annotated.vcf.gz NIST.$var.vcf.gz.tbi
  done
fi

# Prepare the right VCF for SNPs
for var in SNP INDEL
do
  echo "Extracting $var from GATK.vcf.gz"
  (java -Xms1g -Xmx1g -jar $GATK_JAR -T SelectVariants -R $REFERENCE --sample_name NA12878 -selectType $var -V GATK.vcf.gz --excludeFiltered -o CEUTrio.PASS.$var.hg19.vcf -U LENIENT_VCF_PROCESSING &>$LOGDIR/$var.PASS.log; bgzip -f CEUTrio.PASS.$var.hg19.vcf && tabix -f CEUTrio.PASS.$var.hg19.vcf.gz) &
  (java -Xms1g -Xmx1g -jar $GATK_JAR -T SelectVariants -R $REFERENCE --sample_name NA12878 -selectType $var -V GATK.vcf.gz -o CEUTrio.ALL.$var.hg19.vcf -U LENIENT_VCF_PROCESSING &>$LOGDIR/$var.ALL.log; bgzip -f CEUTrio.ALL.$var.hg19.vcf && tabix -f CEUTrio.ALL.$var.hg19.vcf.gz) &
done
wait

for var in SNP INDEL; do
  (
    echo "Extracting $var missed from NIST"
    vcf-isec -c NIST.$var.vcf.gz CEUTrio.PASS.$var.hg19.vcf.gz | bgzip > NIST.$var.missed.in.pass.vcf.gz; tabix -f NIST.$var.missed.in.pass.vcf.gz

    echo "Extracting $var present in ALL, but not in PASS and in NIST"
    vcf-isec CEUTrio.ALL.$var.hg19.vcf.gz NIST.$var.missed.in.pass.vcf.gz | bgzip > NIST.$var.present.in.all.vcf.gz; tabix -f NIST.$var.present.in.all.vcf.gz

    echo "Extracting $var not present in ALL but in NIST"
    vcf-isec -c NIST.$var.missed.in.pass.vcf.gz NIST.$var.present.in.all.vcf.gz | bgzip > NIST.$var.missed.in.all.vcf.gz; tabix -f NIST.$var.missed.in.all.vcf.gz
  ) &
done
wait

vcf-concat CEUTrio.PASS.INDEL.hg19.vcf.gz NIST.INDEL.present.in.all.vcf.gz NIST.INDEL.missed.in.all.vcf.gz | vcf-sort -c | bgzip > INDEL.vcf.gz && tabix -f INDEL.vcf.gz
vcf-concat CEUTrio.PASS.SNP.hg19.vcf.gz NIST.SNP.present.in.all.vcf.gz NIST.SNP.missed.in.all.vcf.gz | vcf-sort -c | bgzip > SNP.vcf.gz && tabix -f SNP.vcf.gz
vcf-concat INDEL.vcf.gz SNP.vcf.gz | vcf-sort -c | bgzip > variants.vcf.gz && tabix -f variants.vcf.gz
vcf-concat deletions.vcf.gz INDEL.vcf.gz | vcf-sort -c | bgzip > SV.INDEL.vcf.gz && tabix -f SV.INDEL.vcf.gz

for vcf in deletions.vcf.gz INDEL.vcf.gz SNP.vcf.gz; do
  (
    vcfunzipped=`basename $vcf .gz`
    gunzip -c $vcf > $vcfunzipped
  ) &
done
wait

cd $GENOMEDIR
java -jar $VCF2DIPLOID -id NA12878 -chr $REFERENCE -vcf $VCFDIR/deletions.vcf $VCFDIR/INDEL.vcf $VCFDIR/SNP.vcf &>$LOGDIR/vcf2diploid.log

echo "Getting list of chromosomes"
all_list=

# Create the bedfile for the chromosomes specified
if [ -z "$CHR_LIST" ]; then
  CHR_LIST=`awk '{print $1}' $REFERENCE.fai`
fi

for chromosome in $CHR_LIST; do
  awk -v chr=$chromosome '{if ($1 == chr) print $1 "\t0\t" $2}' $REFERENCE.fai
done > $GENOMEDIR/regions.bed

for chr in $CHR_LIST; do
  [ -e "$chr"_NA12878_maternal.fa ] && all_list="$all_list $chr""_NA12878_maternal.fa"
  [ -e "$chr"_NA12878_paternal.fa ] && all_list="$all_list $chr""_NA12878_paternal.fa"
done

cat $all_list > NA12878.fa
$SAMTOOLS faidx NA12878.fa
