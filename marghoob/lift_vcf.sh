#!/bin/bash

myname=`basename $0`

function usage {
  echo "$myname <srcid> <dstid> <inputvcf> <outputvcf>"
  exit 1
}

[ $# -ne 4 ] && usage

src=$1
dst=$2
inputvcf=$3
outputvcf=$4

export JAVA_HOME=$HOME/lake/opt/jdk1.7.0_25/
export PATH=$HOME/vcftools_0.1.11/bin:$JAVA_HOME/bin:$PATH
export GATK_DEV=/mnt/scratch0/marghoob/git/gatk-protected/
export CHAIN_DIR=$HOME/lake/users/marghoob/ftp.broadinstitute.org/Liftover_Chain_Files/
export b36=$HOME/lake/users/marghoob/ftp.broadinstitute.org/bundle/2.5/b36/human_b36_both
export b37=$HOME/lake/users/marghoob/ftp.broadinstitute.org/bundle/2.5/b37/human_g1k_v37
export hg19=$HOME/lake/users/marghoob/GATK-bundle-hg19/ucsc.hg19

CHAIN=$CHAIN_DIR/"$src"to"$dst".chain

$GATK_DEV/public/perl/liftOverVCF.pl -vcf $inputvcf -chain $CHAIN -out $outputvcf -gatk $GATK_DEV -newRef ${!dst} -oldRef ${!src} -tmp $PWD
