#!/bin/bash

#$ -S /bin/bash
#$ -j y

[ $# -ne 1 ] && echo "./script jobdir"

jobdir=$1
export PATH="/net/kodiak/volumes/lake/shared/opt/samtools:$PATH"
merge_bams="/net/kodiak/volumes/river/shared/users/marghoob/git/sandbox/marghoob/merge_bams.sh"

MAPPED_BAMS_DIR=$jobdir/bams/recalibrated UNMAPPED_BAMS_DIR=$jobdir/bams/unmapped OUTPUT_BAM=$jobdir/bams/merged.bam SAMTOOLS="/net/kodiak/volumes/lake/shared/opt/samtools/samtools" $merge_bams
