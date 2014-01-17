#!/bin/bash

set -ex

LAKE=/net/kodiak/volumes/lake/shared/
RIVER=/net/kodiak/volumes/river/shared/

INSTALL_PREFIX=$RIVER/users/marghoob/synthetic_genome/bina

BAM2CFG=$INSTALL_PREFIX/breakdancer/current/bin/bam2cfg.pl

function print_abs_path {
  echo $(cd $(dirname $1); pwd)/$(basename $1)
}

myname=`basename $0`
function usage {
  echo "PINDEL_OPTS= WORKDIR=<workdir> LOGDIR= PINDEL=<pindel> BAM_DIR= UNMAPPED_DIR= CHR_LIST= $myname"
  echo "CHR_LIST, PINDEL_OPTS are optional"
  echo "Atleast one of BAM_DIR and UNMAPPED_DIR needs to be specified"
  exit 1
}

[ -z "$WORKDIR" -o -z "$LOGDIR" -o -z "$PINDEL" ] && usage
[ -z "$BAM_DIR" -a -z "$UNMAPPED_DIR" ] && usage
[ ! -e "$PINDEL" ] && echo "$PINDEL not found" && exit 1

mkdir -pv $WORKDIR $LOGDIR
REFERENCE=$LAKE/users/marghoob/GATK-bundle-hg19/ucsc.hg19.fa

if [ -z "$CHR_LIST" ]; then
  CHR_LIST=`awk '{print $1}' $REFERENCE.fai`
fi

WORKDIR=$(print_abs_path $WORKDIR)
[ -n "$BAM_DIR" ] && BAM_DIR=$(print_abs_path $BAM_DIR)
[ -n "$UNMAPPED_DIR" ] && UNMAPPED_DIR=$(print_abs_path $UNMAPPED_DIR)
PINDEL=$(print_abs_path $PINDEL)

START=$(date +%s)

MAPPED_BAMS=
for chr in $CHR_LIST; do
  [ -e "$BAM_DIR/$chr.bam" ] && MAPPED_BAMS="$MAPPED_BAMS $BAM_DIR/$chr.bam"
done
# Generate the pindel config
PERL5LIB=$INSTALL_PREFIX/breakdancer/current/perl $BAM2CFG $MAPPED_BAMS > $WORKDIR/bam2cfg.cfg 2>$LOGDIR/bam2cfg.log
mean_insert_size=`awk '{ split($9, arr_split, ":"); sum = sum + arr_split[2]; } END { print sum / NR; }' $WORKDIR/bam2cfg.cfg`
mkdir -pv $WORKDIR/bams/unmapped $WORKDIR/bams/mapped
rm -f $WORKDIR/bams/unmapped/* $WORKDIR/bams/mapped/*

# Create the symlinks
if [ -n "$UNMAPPED_DIR" ]; then
  for bam in `ls $UNMAPPED_DIR/*.bam`; do
    (cd $WORKDIR/bams/unmapped && ln -sfv $bam && ln -sfv $UNMAPPED_DIR/`basename $bam .bam`.bai `basename $bam`.bai)
  done
fi

for bam in $MAPPED_BAMS; do
  (cd $WORKDIR/bams/mapped && ln -sfv $bam && ln -sfv $BAM_DIR/`basename $bam .bam`.bai `basename $bam`.bai)
done

UNMAPPED_BAMS=
[ -n "$UNMAPPED_DIR" ] && UNMAPPED_BAMS=`ls $WORKDIR/bams/unmapped/*.bam`
for chr in $CHR_LIST; do
  [ ! -s "$BAM_DIR/$chr.bam" ] && continue
  echo "$WORKDIR/bams/mapped/$chr.bam $mean_insert_size pindel_sample" > $WORKDIR/$chr.cfg
  for bam in $UNMAPPED_BAMS; do
    bamid=`basename $bam .bam`
    echo "$bam $mean_insert_size pindel_sample" >> $WORKDIR/$chr.cfg
    #echo "$bam 0 pindel_sample" >> $WORKDIR/$chr.cfg
  done
  
  BREAKDANCER_OPT=
  [ -n "$BREAKDANCER_DIR" ] && BREAKDANCER_OPT="-b $BREAKDANCER_DIR/$chr.out"
  echo "$PINDEL -f $REFERENCE -c $chr -i $WORKDIR/$chr.cfg -o $WORKDIR/$chr.out $BREAKDANCER_OPT $PINDEL_OPTS &>$LOGDIR/$chr.log"
done | xargs -I CMD --max-procs=12 bash -c CMD

END=$(date +%s)
DIFF=$(( $END - $START ))
echo "It took $DIFF seconds" > $LOGDIR/time.log

