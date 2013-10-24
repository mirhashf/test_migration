#!/bin/bash

set -ex

export JAVA_HOME=~/lake/opt/jdk1.7.0_25/
export PATH=$HOME/vcftools_0.1.11/bin:$JAVA_HOME/bin:$PATH:/usr/lib/bina/samtools/current/bin/
CHR_LIST="chr22"

GATK_JAR=$HOME/lake/opt/gatk-2.7-2-g6bda569/GenomeAnalysisTK.jar
BWARUNNER_DIR=$PWD/bina/bwarunner/current/bin
BWARUNNER=$BWARUNNER_DIR/bwarunner.py
BWA=$PWD/bina/bwa/current/bin/bwa
SAMTOOLS=$PWD/bina/samtools/current/bin/samtools
BGZIP=$PWD/bina/tabix/current/bin/bgzip
TABIX=$PWD/bina/tabix/current/bin/tabix
BWAINDEX=$PWD/indexes/bwaindex
HEADER=$PWD/header.sam

myname=`basename $0`

function usage {
  echo "FASTQ1=<reads1> FASTQ2=<reads2> WORKDIR=<workdir> $myname"
  exit 1
}

[ -z "$FASTQ1" ] && usage
[ -z "$FASTQ2" ] && usage
[ -z "$WORKDIR" ] && usage

WORKDIR=$PWD/$WORKDIR
LOGDIR=$WORKDIR/logs

mkdir -pv $WORKDIR
mkdir -pv $LOGDIR

dbsnp=~/lake/users/marghoob/GATK-bundle-hg19/dbsnp_137.hg19.vcf
reference=~/lake/users/marghoob/GATK-bundle-hg19/ucsc.hg19.fa

pushd $PWD
cd $WORKDIR

JAVA_TMP=$WORKDIR/tmp

INTERVALS=
REGIONS=
for chr in $CHR_LIST; do
  awk -v chr=$chr '{if ($1 == chr) { print $1" 0 "$2 }}' $reference.fai > $chr.bed
  INTERVALS="$INTERVALS -L $chr.bed"
  if [ -z "$REGIONS" ]; then
    REGIONS="-r $chr"
  else
    REGIONS="$REGIONS,$chr"
  fi
done

cat $HEADER <($BWARUNNER --bwa $BWA --splitter $BWARUNNER_DIR/fastq_splitter --merger $BWARUNNER_DIR/sam_merger --reference $BWAINDEX --reads1 $FASTQ1 --reads2 $FASTQ2 --bwa_aln_options '-t 4 -q 30' --bwa_samse_options ' -r "@RG	ID:lane0	SM:NA12878	LB:pairedend"' --bwa_sampe_options '-P -r "@RG	ID:lane0	SM:NA12878	LB:pairedend"' --log_directory $LOGDIR/bwarunner_log.0 --working_directory $WORKDIR/bwarunner_work --compression_level 1 --nthreads 8 --splitter_blocksize '65536' --splitter_max_pending '4' 2>$LOGDIR/aligner.0.log) | $SAMTOOLS view -b -1 -S - > bwa.aligned.bam
$SAMTOOLS sort -m 50000000000 bwa.aligned.bam bwa.aligned.sorted
$SAMTOOLS index bwa.aligned.sorted.bam

# Create the bed file

java -Djava.io.tmpdir=$JAVA_TMP -jar $GATK_JAR -T RealignerTargetCreator --input_file bwa.aligned.sorted.bam --out realign.intervals --reference_sequence $reference $INTERVALS 2>&1 1>$LOGDIR/realigner_target_creator.log

java -Djava.io.tmpdir=$JAVA_TMP -jar $GATK_JAR -T IndelRealigner --input_file bwa.aligned.sorted.bam --out realigned.bam --reference_sequence $reference --targetIntervals realign.intervals 2>&1 1>$LOGDIR/indel_realigner.log
$SAMTOOLS index realigned.bam

java -Djava.io.tmpdir=$JAVA_TMP -jar $GATK_JAR -T BaseRecalibrator -cov 'ReadGroupCovariate' -cov 'QualityScoreCovariate' -cov 'CycleCovariate' -cov 'ContextCovariate' --input_file realigned.bam $INTERVALS --reference_sequence $reference -o recal.table -knownSites $dbsnp &>$LOGDIR/base_recalibrator.log

java -Djava.io.tmpdir=$JAVA_TMP -jar $GATK_JAR -T PrintReads --input_file realigned.bam --out recalibrated.bam --reference_sequence $reference -BQSR recal.table -baq RECALCULATE $INTERVALS &>$LOGDIR/print_reads.log
$SAMTOOLS index recalibrated.bam

java -Djava.io.tmpdir=$JAVA_TMP -jar $GATK_JAR -T UnifiedGenotyper --input_file recalibrated.bam --out genotyped.vcf --reference_sequence $reference --dbsnp $dbsnp $INTERVALS -stand_emit_conf 0.1 -baq CALCULATE_AS_NECESSARY --genotype_likelihoods_model BOTH -A AlleleBalance -A Coverage -A MappingQualityZero -nt 8 &>$LOGDIR/unified_genotyper.log

$BGZIP -f genotyped.vcf
$TABIX -f genotyped.vcf.gz

for vartype in SNP INDEL; do
  (java -Djava.io.tmpdir=$JAVA_TMP -jar $GATK_JAR -T SelectVariants -V genotyped.vcf.gz -o genotyped.$vartype.vcf -selectType $vartype -R $reference &>$LOGDIR/SelectVariants.$vartype.log; $BGZIP -f genotyped.$vartype.vcf && $TABIX -f genotyped.$vartype.vcf.gz) &
done
wait

for vartype in SNP INDEL; do
  vcf-compare $REGIONS genotyped.$vartype.vcf.gz /home/marghoob/lake/users/marghoob/NIST/NISThighConf.$vartype.hg19.vcf.gz > NIST.sensitivity.$vartype.txt
done

popd
