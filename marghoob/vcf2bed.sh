#!/bin/bash

function usage {
    echo "convert2bed.sh <vcf.gz> <out.bed>"
    exit 1
}

vcfgz=$1
outfile=$2

[ -z "$vcfgz" ] && usage
[ -z "$outfile" ] && usage

gunzip -c $vcfgz|awk '{OFS="\t"; if (!/^#/){print $1,$2-1,$2,$4"/"$5,"+"}}' > $outfile
