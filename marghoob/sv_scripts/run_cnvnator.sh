#!/bin/bash

set -e

LAKE=/net/kodiak/volumes/lake/shared
RIVER=/net/kodiak/volumes/river/shared
export JAVA_HOME=$LAKE/opt/jdk1.7.0_25/
export PATH=$HOME/vcftools_0.1.11/bin:$JAVA_HOME/bin:$PATH:/usr/lib/bina/samtools/current/bin/:$PWD/bina/bwa/current/bin/
CHR_LIST=

LAKE=$HOME/lake

GATK_JAR=$LAKE/opt/gatk-2.7-2-g6bda569/GenomeAnalysisTK.jar
INSTALL_PREFIX=$RIVER/synthetic_genome/bina
export CONTIGS=$RIVER/synthetic_genome/indexes/b37/contigs

myname=`basename $0`

function usage {
  echo "WORKDIR=<workdir> BAMS_DIR= LOGDIR= $myname"
  exit 1
}

function print_abs_path {
  echo $(cd $(dirname $1); pwd)/$(basename $1)
}

[ -z "$WORKDIR" -o -z "$BAMS_DIR" -o -z "$LOGDIR" ] && usage

mkdir -pv $WORKDIR $LOGDIR $WORKDIR/bams

WORKDIR=$(print_abs_path $WORKDIR)
LOGDIR=$(print_abs_path $LOGDIR)
BAMS_DIR=$(print_abs_path $BAMS_DIR)

REFERENCE=~/lake/users/marghoob/GATK-bundle-hg19/ucsc.hg19.fa

START=$(date +%s)

if [ -z "$CHR_LIST" ]; then
  CHR_LIST=`awk '{print $1}' $REFERENCE.fai`
fi

for chr in $CHR_LIST; do
  [ ! -e "$BAMS_DIR/$chr.bam" ] && continue
  (cd $WORKDIR/bams; ln -sf $BAMS_DIR/$chr.bam; ln -sf $BAMS_DIR/$chr.bai $chr.bam.bai)
done

for chr in $CHR_LIST; do
  echo "INSTALL_PREFIX=$INSTALL_PREFIX OUTPUT_PREFIX=$WORKDIR LOG_PREFIX=$LOGDIR SCRATCH_PREFIX=$WORKDIR  CHROMOSOME=$chr VERSION=current CONTIGS=$CONTIGS BIN_SIZE=100 $INSTALL_PREFIX/cnvnator/current/bin/cnvnator-chromosome-runner $WORKDIR/bams/$chr.bam"
done | xargs -I CMD --max-procs=24 bash -c CMD 

END=$(date +%s)
DIFF=$(( $END - $START ))
echo "It took $DIFF seconds" > $LOGDIR/time.log
