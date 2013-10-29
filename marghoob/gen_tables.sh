#!/bin/bash

set -ex

indir=$1
label1=$2
label2=$3

declare -A labels
labels["common"]="$label1 and $label2"
labels["1"]="$label1 only"
labels["2"]="$label2 only"

# Get counts for the SNP tables
for filter in ALL PASS; do
  for subset in SNP INDEL; do
    for subsubset in common 1 2; do
      basedir="$indir/$filter/$subset"
      total=`cat $basedir/$subsubset.stats.counts | grep -v "^#" | head -n 1 | awk '{print $1}'`
      titv=`cat $basedir/$subsubset.stats.tstv | grep -v "^#" | head -n 1 | awk '{print $3}'`
      total_known=`cat $basedir/$subsubset/known.counts | grep -v "^#" | head -n 1 | awk '{print $1}'`
      known_titv=`cat $basedir/$subsubset/known.tstv | grep -v "^#" | head -n 1 | awk '{print $3}'`
      novel_titv=`cat $basedir/$subsubset/novel.tstv | grep -v "^#" | head -n 1 | awk '{print $3}'`
      known_ratio=`echo "scale=2; 100.0*$total_known/$total"|bc -l`
      novel_ratio=`echo "scale=2; 100.0 - $known_ratio"|bc -l`

      echo "$filter	${labels[$subsubset]}"
      echo "Sites	$subset	Ti/Tv"
      echo "All	$total	$titv"
      echo "Known	$known_ratio	$known_titv"
      echo "Novel	$novel_ratio	$novel_titv"
      echo ""
    done
  done
done
