#!/bin/bash

export JAVA_HOME=$HOME/lake/opt/jdk1.7.0_25/
export PATH=$HOME/vcftools_0.1.11/bin:$JAVA_HOME/bin:$PATH
export GATK_DEV=/mnt/scratch0/marghoob/git/gatk-protected/
export CHAIN_DIR=$HOME/lake/users/marghoob/ftp.broadinstitute.org/Liftover_Chain_Files/
export B36=$HOME/lake/users/marghoob/ftp.broadinstitute.org/bundle/2.5/b36/human_b36_both
export B37=$HOME/lake/users/marghoob/ftp.broadinstitute.org/bundle/2.5/b37/human_g1k_v37
export HG19=$HOME/lake/users/marghoob/GATK-bundle-hg19/ucsc.hg19

input=$1
output=$2

$GATK_DEV/public/perl/liftOverVCF.pl -vcf $input -chain $CHAIN_DIR/b37tohg19.chain -out $output -gatk $GATK_DEV -newRef $HG19 -oldRef $B37 -tmp $PWD
