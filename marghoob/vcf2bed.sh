#!/bin/bash

function usage {
    echo "convert2bed.sh <vcf.gz> <out.bed>"
    exit 1
}

[ $# -ne 2 ] && usage

vcfgz=$1
outfile=$2

gunzip -c $vcfgz|awk '{OFS="\t"; if (!/^#/){print $1,$2-1,$2,$4"/"$5,"+"}}' > $outfile
