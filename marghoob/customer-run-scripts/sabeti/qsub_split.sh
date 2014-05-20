#!/bin/bash

#$ -S /bin/bash
#$ -j y

[ $# -ne 1 ] && echo "./script sampledir"

export PATH="/net/kodiak/volumes/lake/shared/opt/samtools:$PATH"
split_bam=/net/kodiak/volumes/river/shared/users/marghoob/git/sandbox/marghoob/split_bam/split_bam
samples=$1
for sample in $samples; do
  echo "--------------------------------------------------------" >&2
  echo $sample >&2
  bams=`ls $sample/*.bam`
  #rgs=
  for bam in $bams; do
    [[ $bam = *.rg.bam ]] && echo "$bam will be skipped" >&2 && continue
    pathname=`dirname $bam`
    prefix=`basename $bam .bam`
    LD_LIBRARY_PATH=/net/kodiak/volumes/river/shared/users/marghoob/git/htslib $split_bam $bam $pathname/$prefix
  done
  for bam in `ls $sample/*.rg.bam`; do
    samtools index $bam &
  done
  #rgs=`echo $rgs | xargs -n1 | sort -u | xargs`
done
wait
