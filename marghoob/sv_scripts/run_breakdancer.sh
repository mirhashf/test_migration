#!/bin/bash

set -ex

export JAVA_HOME=~/lake/opt/jdk1.7.0_25/
export PATH=$HOME/vcftools_0.1.11/bin:$JAVA_HOME/bin:$PATH:/usr/lib/bina/samtools/current/bin/:$PWD/bina/bwa/current/bin/

LAKE=$HOME/lake

BAM2CFG=$PWD/bina/breakdancer/current/bin/bam2cfg.pl
BREAKDANCER=$PWD/bina/breakdancer/current/bin/breakdancer_max
INSTALL_PREFIX=$PWD/bina
SAMTOOLS=$LAKE/opt/samtools/samtools

myname=`basename $0`

function usage {
  echo "BAMS_DIR= WORKDIR=<workdir> LOGDIR= OPTS= $myname"
  exit 1
}

function print_abs_path {
  echo $(cd $(dirname $1); pwd)/$(basename $1)
}

[ -z "$WORKDIR" -o -z "$BAMS_DIR" -o -z "$LOGDIR" ] && usage

REFERENCE=/net/kodiak/volumes/lake/shared/users/marghoob/GATK-bundle-hg19/ucsc.hg19.fa

WORKDIR=$(print_abs_path $WORKDIR)
LOGDIR=$(print_abs_path $LOGDIR)
BAMS_DIR=$(print_abs_path $BAMS_DIR)

START=$(date +%s)

if [ -z "$CHR_LIST" ]; then
  CHR_LIST=`awk '{print $1}' $REFERENCE.fai`
fi

mkdir -pv $WORKDIR $LOGDIR $WORKDIR/bams

for chr in $CHR_LIST; do
  [ ! -e "$BAMS_DIR/$chr.bam" ] && continue
  (cd $WORKDIR/bams; ln -sf $BAMS_DIR/$chr.bam; ln -sf $BAMS_DIR/$chr.bai $chr.bam.bai) > /dev/null
  echo "PERL5LIB=$INSTALL_PREFIX/breakdancer/current/perl $BAM2CFG $WORKDIR/bams/$chr.bam 2> $LOGDIR/$chr.bam2cfg.log > $WORKDIR/$chr.cfg"
done | xargs -I CMD --max-procs=24 bash -c CMD

for chr in $CHR_LIST; do
  [ ! -e "$BAMS_DIR/$chr.bam" ] && continue
  echo "$BREAKDANCER $OPTS -o $chr $WORKDIR/$chr.cfg 2>$LOGDIR/$chr.log > $WORKDIR/$chr.out"
done | xargs -I CMD --max-procs=24 bash -c CMD

END=$(date +%s)
DIFF=$(( $END - $START ))
echo "It took $DIFF seconds" > $LOGDIR/time.log

