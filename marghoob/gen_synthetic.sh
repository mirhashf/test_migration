#!/bin/bash

set -ex

myname=`basename $0`
function usage {
  echo "START_E= END_E= $myname"
  exit 1
}

NIST=$HOME/lake/users/marghoob/NIST/NISThighConf
GATKVCF=~/lake/users/marghoob/GATK-bundle-hg19/CEUTrio.HiSeq.WGS.b37.bestPractices.phased.hg19.vcf.gz
DELETIONS_VCF=~/lake/users/marghoob/homes.gersteinlab.org/people/aabyzov/outbox/private/NA12878_variants/NA12878.2010_06.and.fosmid.deletions.phased.vcf

#CHR_LIST="chr22"
#START_E=0.0001
#END_E=0.001

[ -z "$START_E" ] && usage
[ -z "$END_E" ] && usage
#[ -z "$CHR_LIST" ] && usage

SCRATCHDIR=$HOME/scratch/dwgsim
WORKDIR=$PWD/work
LOGDIR=$WORKDIR/log
mkdir -p $WORKDIR
mkdir -pv $SCRATCHDIR
mkdir -pv $LOGDIR

# Tool vars
LIFTVCF=~/git/sandbox/marghoob/lift_vcf.sh
GATK_JAR=/home/marghoob/lake/opt/CancerAnalysisPackage-2013.2-18-g8207e53/GenomeAnalysisTK.jar
VCF2DIPLOID=~/lake/users/marghoob/vcf2diploid/vcf2diploid.jar
IGVTOOLS=~/lake/users/marghoob/IGVTools/igvtools.jar
DWGSIM=~/lake/opt/dwgsim/dwgsim
BEDTOOLS_DIR=~/lake/opt/bedtools-2.17.0/bin/
SAMTOOLS=/usr/lib/bina/samtools/current/bin/samtools

REFERENCE=~/lake/users/marghoob/GATK-bundle-hg19/ucsc.hg19.fa

cd $WORKDIR

# Create the bedfile for the chromosomes specified
if [ -z "$CHR_LIST" ]; then
  CHR_LIST=`awk '{print $1}' $REFERENCE.fai`
fi

if (( 0 )); then

# Prepare the right VCF for deletions
echo "Lifting the deletion SVs from b36 to b37 after removing records with incorrect filter and removing the END field in INFO column"
awk '/^#/ { str = $0; sub(/-1/, "1", str); print str} !/^#/ { str = $0; sub(/END=[0-9]+;/, "", str); if ($7 != "0") print str}' $DELETIONS_VCF > deletions.b36.vcf
$LIFTVCF b36 b37 deletions.b36.vcf deletions.b37.vcf &>$LOGDIR/deletions.b36.to.b37.log
echo "Lifting the deletion SVs from b37 to hg19"
$LIFTVCF b37 hg19 deletions.b37.vcf deletions.hg19.vcf &>$LOGDIR/deletions.b37.to.hg19.log

# Prepare the right VCF for SNPs
for var in SNP INDEL
do
  echo "Extracting $var from $GATKVCF"
  (java -Xms1g -Xmx1g -jar $GATK_JAR -T SelectVariants -R $REFERENCE --sample_name NA12878 -selectType $var -V $GATKVCF --excludeFiltered -o CEUTrio.PASS.$var.hg19.vcf -U LENIENT_VCF_PROCESSING &>$LOGDIR/$var.PASS.log; bgzip -f CEUTrio.PASS.$var.hg19.vcf && tabix -f CEUTrio.PASS.$var.hg19.vcf.gz) &
  (java -Xms1g -Xmx1g -jar $GATK_JAR -T SelectVariants -R $REFERENCE --sample_name NA12878 -selectType $var -V $GATKVCF -o CEUTrio.ALL.$var.hg19.vcf -U LENIENT_VCF_PROCESSING &>$LOGDIR/$var.ALL.log; bgzip -f CEUTrio.ALL.$var.hg19.vcf && tabix -f CEUTrio.ALL.$var.hg19.vcf.gz) &
done
wait

for var in SNP INDEL; do
  (echo "Extracting $var missed from NIST"
  vcf-isec -c $NIST.$var.hg19.annotated.vcf.gz CEUTrio.PASS.$var.hg19.vcf.gz | bgzip > NIST.$var.missed.in.pass.vcf.gz; tabix -f NIST.$var.missed.in.pass.vcf.gz

  echo "Extracting $var present in ALL, but not in PASS and in NIST"
  vcf-isec CEUTrio.ALL.$var.hg19.vcf.gz NIST.$var.missed.in.pass.vcf.gz | bgzip > NIST.$var.present.in.all.vcf.gz; tabix -f NIST.$var.present.in.all.vcf.gz

  echo "Extracting $var not present in ALL but in NIST"
  vcf-isec -c NIST.$var.missed.in.pass.vcf.gz NIST.$var.present.in.all.vcf.gz | bgzip > NIST.$var.missed.in.all.vcf.gz; tabix -f NIST.$var.missed.in.all.vcf.gz) &
done
wait

for vcf in *.vcf.gz; do
  (gunzip $vcf && java -Xms1g -Xmx1g -jar $IGVTOOLS index $vcf) &
done
wait
fi

# Generate the bedfiles from the VCFs
if (( 0 )); then
awk '!/^#/ {
             split($8, info_split, ";");
             svlen = 0;
             for (i in info_split) {
               if (match(info_split[i], "^SVLEN=")) {
                 split(info_split[i], field_split, "=");
                 svlen = -field_split[2];
               }
             };
             print $1"\t"($2-1)"\t"($2+svlen-1);
           }' deletions.hg19.vcf > deletions.hg19.bed

for SNP_vcf in CEUTrio.PASS.SNP.hg19.vcf NIST.SNP.present.in.all.vcf NIST.SNP.missed.in.all.vcf; do
  prefix=`basename $SNP_vcf .vcf`
  awk '!/^#/ {print $1"\t"($2-1)"\t"$2}' $SNP_vcf > $prefix.bed
done
fi

# Build the genome
mkdir -pv vcf2diploid
cd vcf2diploid
#java -jar $VCF2DIPLOID -id NA12878 -chr $REFERENCE -vcf ../deletions.hg19.vcf ../CEUTrio.PASS.INDEL.hg19.vcf ../NIST.INDEL.present.in.all.vcf ../NIST.INDEL.missed.in.all.vcf ../CEUTrio.PASS.SNP.hg19.vcf ../NIST.SNP.present.in.all.vcf ../NIST.SNP.missed.in.all.vcf &>$LOGDIR/vcf2diploid.log

echo "Getting list of chromosomes"
all_list=

for chr in $CHR_LIST; do
  [ -e "$chr"_NA12878_maternal.fa ] && all_list="$all_list $chr""_NA12878_maternal.fa"
  [ -e "$chr"_NA12878_paternal.fa ] && all_list="$all_list $chr""_NA12878_paternal.fa"
done

cat $all_list > NA12878.fa
$SAMTOOLS faidx NA12878.fa

cd ..

mkdir -pv fastq
for lane in 0 1 2; do
  $DWGSIM -e $START_E,$END_E -E $START_E,$END_E -r 0 -F 0 -d 330 -s 70 -C 6 -1 100 -2 100 -z $lane vcf2diploid/NA12878.fa fastq/simulated.lane$lane

  (gzip -1 -f fastq/simulated.lane$lane.bwa.read1.fastq)&
  (gzip -1 -f fastq/simulated.lane$lane.bwa.read2.fastq)&
  wait
done
