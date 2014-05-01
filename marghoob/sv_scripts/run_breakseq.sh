#!/bin/bash

set -e

export PATH=$HOME/vcftools_0.1.11/bin:$JAVA_HOME/bin:$PATH:/usr/lib/bina/samtools/current/bin/:$PWD/bina/bwa/current/bin/

LAKE=$HOME/lake

BWA=$PWD/bina/bwa/current/bin/bwa
SAMTOOLS=$LAKE/opt/samtools/samtools

myname=`basename $0`

function usage {
  echo "BREAKSEQ= BPLIB= UNMAPPED_DIR= UNMAPPED_BAMS= MAPPED_DIR= MAPPED_BAMS= BWA= SAMTOOLS= WORKDIR= LOGDIR= $myname"
  echo "Need either one of UNMAPPED_DIR or UNMAPPED_BAMS to be specified"
  exit 1
}

[ -z "$WORKDIR" -o -z "$BPLIB" -o -z "$BWA" -o -z "$SAMTOOLS" ] && usage
[ -z "$UNMAPPED_DIR" -a -z "$UNMAPPED_BAMS" -a -z "$MAPPED_BAMS" ] && usage
[ -z "$MIN_SOFT_CLIP" ] && echo "MIN_SOFT_CLIP not specified" && exit 1
[ -z "$MIN_SOFT_CLIP_MAPQ" ] && echo "MIN_SOFT_CLIP_MAPQ not specified" && exit 1
[ -z "$MIN_SOFT_CLIP_MATE_MAPQ" ] && echo "MIN_SOFT_CLIP_MATE_MAPQ not specified" && exit 1
[ -z "$BAD_MAP_MAX_SOFT_CLIP" ] && echo "BAD_MAP_MAX_SOFT_CLIP not specified" && exit 1
[ -z "$BAD_MAP_MIN_MAPQ" ] && echo "BAD_MAP_MIN_MAPQ not specified" && exit 1
[ -z "$BAD_MAP_MIN_NM" ] && echo "BAD_MAP_MIN_NM not specified" && exit 1
[ -z "$BAD_MAP_MIN_MATE_MAPQ" ] && echo "BAD_MAP_MIN_MATE_MAPQ not specified" && exit 1

mkdir -pv $WORKDIR $LOGDIR

[ -n "$UNMAPPED_DIR" ] && UNMAPPED_BAMS=`ls $UNMAPPED_DIR/*.bam`
[ -n "$MAPPED_DIR" ] && MAPPED_BAMS=`ls $MAPPED_DIR/*.bam`
BAMS="$MAPPED_BAMS $UNMAPPED_BAMS"

START=$(date +%s)

MIN_SOFT_CLIP_MATE_MAPQ=$MIN_SOFT_CLIP_MATE_MAPQ BAD_MAP_MIN_MATE_MAPQ=$BAD_MAP_MIN_MATE_MAPQ MIN_SOFT_CLIP_MAPQ=$MIN_SOFT_CLIP_MAPQ BAD_MAP_MAX_SOFT_CLIP=$BAD_MAP_MAX_SOFT_CLIP BAD_MAP_MIN_MAPQ=$BAD_MAP_MIN_MAPQ BAD_MAP_MIN_NM=$BAD_MAP_MIN_NM MIN_SOFT_CLIP="$MIN_SOFT_CLIP" BPLIB="$BPLIB" SAMTOOLS="$SAMTOOLS" BWA="$BWA" LD_LIBRARY_PATH=$BREAKSEQ/lib $BREAKSEQ/bin/breakseq $WORKDIR/breakseq_out.gff $BAMS &>$LOGDIR/breakseq.log

END=$(date +%s)
DIFF=$(( $END - $START ))
echo "It took $DIFF seconds" > $LOGDIR/time.log
