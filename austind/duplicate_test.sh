#!/bin/bash

usage ( ) {
    echo "usage: $0 in.sam in.bed"
    exit 1
}

if [ -z "$1" ]
then
    usage
fi

if [ -z "$2" ]
then
    usage
fi

INPUT_SAM="$1"
INPUT_BED="$2"

ROOT=$(pwd)/$(dirname $0)

echo "starting samtools mapped filter..."
samtools view -h -F 0x4 -f 0x2 -S $INPUT_SAM > $INPUT_SAM.mapped.sam 2> /dev/null

SAMSORTER_BIN="$ROOT/../build/samsorter"

echo "starting samsorter..."
mkdir -p spill depth_of_coverage ui_coverage
$SAMSORTER_BIN --sam_filename $INPUT_SAM.mapped.sam --bed_filename $INPUT_BED --max_mem_gb 8 --remove_duplicates 2> /dev/null
rm -rf spill depth_of_coverage ui_coverage

PICARD_JAR="$ROOT/../../third-party/picard-1.64.jar"
SAMTOOLS_JAR="$ROOT/../../third-party/samtools-1.64.jar"

JP="java -cp $PICARD_JAR:$SAMTOOLS_JAR"

O=$(basename $INPUT_SAM).mapped.csorted

echo "starting picard sort..."
$JP net.sf.picard.sam.SortSam I=$INPUT_SAM.mapped.sam O=$O.bam SORT_ORDER=coordinate 2> /dev/null

echo "starting picard duplicate marking..."
$JP net.sf.picard.sam.MarkDuplicates I=$O.bam O=$O.marked.bam M=$O.marked.m 2> /dev/null

echo "starting samtools duplicate filter..."
samtools view -h -f 0x400 -b $O.marked.bam > $O.duplicates.bam

ALIGNSTATS_JAR="$ROOT/../../tools/alignstats/target/alignstats.jar"

ALIGNSTATS_OPTS="--bins --report --lenient-validation-stringency --names samsorter,picard"

echo "starting alignstats..."
java -jar $ALIGNSTATS_JAR out_0_duplicates.bam $O.duplicates.bam $ALIGNSTATS_OPTS 2> /dev/null
