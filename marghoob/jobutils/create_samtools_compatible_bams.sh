#!/bin/bash

jobdir=$1
symlinkdir=$2

[ -z "$jobdir" -o -z "$symlinkdir" ] && echo "Usage: ./script.sh <jobdir> <symlinkdir>" && exit 1

jobdir=`readlink -f $jobdir`

mkdir -pv $symlinkdir
cd $symlinkdir
for bam in `ls $jobdir/bams/recalibrated/*.bam $jobdir/bams/unmapped/*.bam`; do
  prefix=`basename $bam .bam`
  ln -sv $bam
  ln -sv "${bam%.bam}.bai" $prefix.bam.bai
done
