#!/bin/bash

LAKE=/net/kodiak/volumes/lake/shared/opt
export PATH=$LAKE/samtools:$PATH

function usage {
  echo "./script <bam>"
  exit 1
}

[ $# -lt 1 ] && usage

samtools view -H $1 |grep "^@RG"|awk '{print substr($2,4)}'
