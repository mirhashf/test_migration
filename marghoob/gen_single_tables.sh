#!/bin/bash

indir=$1
label=$2

echo "Stats for $label"
# Get counts for the SNP tables
for subset in SNP INDEL
do

total=`cat $indir/$subset.counts | grep -v "^#" | head -n 1 | awk '{print $1}'`
titv=`cat $indir/$subset.tstv | grep -v "^#" | head -n 1 | awk '{print $3}'`
total_known=`cat $indir/$subset/known.counts | grep -v "^#" | head -n 1 | awk '{print $1}'`
total_novel=`cat $indir/$subset/novel.counts | grep -v "^#" | head -n 1 | awk '{print $1}'`
known_titv=`cat $indir/$subset/known.tstv | grep -v "^#" | head -n 1 | awk '{print $3}'`
novel_titv=`cat $indir/$subset/novel.tstv | grep -v "^#" | head -n 1 | awk '{print $3}'`
known_ratio=`echo "scale=2; 100.0*$total_known/$total"|bc -l`
novel_ratio=`echo "scale=2; 100.0 - $known_ratio"|bc -l`

echo "Sites	$subset	Ti/Tv"
echo "All	$total	$titv"
echo "Known	$total_known ($known_ratio %)	$known_titv"
echo "Novel	$total_novel ($novel_ratio %)	$novel_titv"
echo ""

done
