#!/bin/bash

set -e

myname=`basename $0`
function usage {
  echo "START_E= END_E= NLANES= TOTAL_COVERAGE= CHR_LIST= $myname"
  echo "START_E=Error rate on the first base."
  echo "END_E=Error rate on the last base."
  echo "NLANES=Number of lanes. NLANES pairs of FASTQs will be generated."
  echo "TOTAL_COVERAGE=Total coverage across all lanes."
  echo "CHR_LIST is optional. If it is empty, all chromosomes are simulated."
  exit 1
}

LAKE=/net/kodiak/volumes/lake/shared
RIVER=/net/kodial/volumes/river/shared

NIST=$LAKE/users/marghoob/NIST/NISThighConf
GATKVCF=$LAKE/users/marghoob/GATK-bundle-hg19/CEUTrio.HiSeq.WGS.b37.bestPractices.phased.hg19.vcf.gz
DELETIONS_VCF=$LAKE/users/marghoob/homes.gersteinlab.org/people/aabyzov/outbox/private/NA12878_variants/NA12878.2010_06.and.fosmid.deletions.phased.vcf

[ -z "$START_E" -o -z "$END_E" -o -z "$NLANES" -o -z "$TOTAL_COVERAGE" ] && usage

SCRATCHDIR=$HOME/scratch/dwgsim
WORKDIR=$PWD/work
LOGDIR=$WORKDIR/log
mkdir -pv $WORKDIR $SCRATCHDIR $LOGDIR

# Tool vars
DIR="$( cd "$( dirname "$0" )" && pwd )"
LIFTVCF=$DIR/lift_vcf.sh
GATK_JAR=$LAKE/opt/CancerAnalysisPackage-2013.2-18-g8207e53/GenomeAnalysisTK.jar
VCF2DIPLOID=$LAKE/users/marghoob/vcf2diploid/vcf2diploid.jar
IGVTOOLS=$LAKE/users/marghoob/IGVTools/igvtools.jar
DWGSIM=$LAKE/opt/dwgsim/dwgsim
BEDTOOLS_DIR=$LAKE/opt/bedtools-2.17.0/bin/
SAMTOOLS=$LAKE/opt/samtools/samtools
export PATH=$LAKE/opt/tabix-0.2.6:$LAKE/opt/vcftools_0.1.11/bin:$PATH
export PERL5LIB=$LAKE/opt/vcftools_0.1.11/perl

REFERENCE=$LAKE/users/marghoob/GATK-bundle-hg19/ucsc.hg19.fa

cd $WORKDIR

# Create the bedfile for the chromosomes specified
if [ -z "$CHR_LIST" ]; then
  CHR_LIST=`awk '{print $1}' $REFERENCE.fai`
fi

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
  (
    echo "Extracting $var missed from NIST"
    vcf-isec -c $NIST.$var.hg19.annotated.vcf.gz CEUTrio.PASS.$var.hg19.vcf.gz | bgzip > NIST.$var.missed.in.pass.vcf.gz; tabix -f NIST.$var.missed.in.pass.vcf.gz

    echo "Extracting $var present in ALL, but not in PASS and in NIST"
    vcf-isec CEUTrio.ALL.$var.hg19.vcf.gz NIST.$var.missed.in.pass.vcf.gz | bgzip > NIST.$var.present.in.all.vcf.gz; tabix -f NIST.$var.present.in.all.vcf.gz

    echo "Extracting $var not present in ALL but in NIST"
    vcf-isec -c NIST.$var.missed.in.pass.vcf.gz NIST.$var.present.in.all.vcf.gz | bgzip > NIST.$var.missed.in.all.vcf.gz; tabix -f NIST.$var.missed.in.all.vcf.gz
  ) &
done
wait

for vcf in *.vcf.gz; do
  (
    vcfunzipped=`basename $vcf .gz`
    gunzip -c $vcf > $vcfunzipped && java -Xms1g -Xmx1g -jar $IGVTOOLS index $vcfunzipped
  ) &
done
wait

# Generate the bedfiles from the VCFs
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


# Build the genome
mkdir -pv vcf2diploid
cd vcf2diploid
java -jar $VCF2DIPLOID -id NA12878 -chr $REFERENCE -vcf ../deletions.hg19.vcf ../CEUTrio.PASS.INDEL.hg19.vcf ../NIST.INDEL.present.in.all.vcf ../NIST.INDEL.missed.in.all.vcf ../CEUTrio.PASS.SNP.hg19.vcf ../NIST.SNP.present.in.all.vcf ../NIST.SNP.missed.in.all.vcf &>$LOGDIR/vcf2diploid.log

echo "Getting list of chromosomes"
all_list=

for chr in $CHR_LIST; do
  [ -e "$chr"_NA12878_maternal.fa ] && all_list="$all_list $chr""_NA12878_maternal.fa"
  [ -e "$chr"_NA12878_paternal.fa ] && all_list="$all_list $chr""_NA12878_paternal.fa"
done

cat $all_list > NA12878.fa
$SAMTOOLS faidx NA12878.fa

cd ..

FASTQDIR=fastq.$START_E"_"$END_E"_"$NLANES"_"$TOTAL_COVERAGE
mkdir -pv $FASTQDIR
COVERAGE_PER_LANE=`echo "scale = 2; $TOTAL_COVERAGE / 2.0 / $NLANES"|bc -l`
let NLANES_MINUS=NLANES-1

rm -f $FASTQDIR/README
for lane in `seq 0 $NLANES_MINUS`; do
  echo "$DWGSIM -e $START_E,$END_E -E $START_E,$END_E -y 0.02 -r 0 -F 0 -d 330 -s 70 -C $COVERAGE_PER_LANE -1 100 -2 100 -z $lane vcf2diploid/NA12878.fa $FASTQDIR/simulated.lane$lane" >> $FASTQDIR/README
  $DWGSIM -e $START_E,$END_E -E $START_E,$END_E -y 0.02 -r 0 -F 0 -d 330 -s 70 -C $COVERAGE_PER_LANE -1 100 -2 100 -z $lane vcf2diploid/NA12878.fa $FASTQDIR/simulated.lane$lane
  gzip -v -1 -f $FASTQDIR/simulated.lane$lane.bwa.read1.fastq &>$LOGDIR/$lane.read1.gzip.log &
  gzip -v -1 -f $FASTQDIR/simulated.lane$lane.bwa.read2.fastq &>$LOGDIR/$lane.read2.gzip.log &
  wait
done
