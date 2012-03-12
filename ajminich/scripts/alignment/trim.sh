#!/bin/bash -eu

TRIMMER="/home/ajminich/programs/sandbox/ajminich/scripts/trimBWAstyle.pl"

# Trims reads.

if [[ $# -lt 2 ]]; then
    echo "Read Trimmer"
    echo "Usage: trimmer.sh <reads_prefix> <quality>"
    echo "Trims reads in the BWA style and writes the FASTQ files to <reads_prefix>_trimmed_#.fq."
    echo "  <reads_prefix>  - the read prefix (<prefix>_1.fq and <prefix>_2.fq)"
    echo "  <quality>       - the quality score to trim to"
else
    reads=${1}
    quality=${2}
    
    cat ${reads}_1.fq | perl ${TRIMMER} -q ${quality} > ${reads}\_trimmed\_1.fq
    cat ${reads}_2.fq | perl ${TRIMMER} -q ${quality} > ${reads}\_trimmed\_2.fq
    
fi
