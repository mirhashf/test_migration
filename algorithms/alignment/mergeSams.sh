#!/bin/bash -eu

# Merges SAM files in the current directory into a single BAM file.

if [[ $# -lt 1 ]]; then
    echo "SAM Merger"
    echo "Usage: mergeSams.sh <finalBamFile> [<samPrefix>]"
    echo "Merges SAM files into a single BAM file."
    echo "  <finalBamFile>   - the final BAM file to output"
    echo "  <samPrefix>      - the prefix for the SAM files. If no prefix"
    echo "                     is provided, all SAM files will be merged."
    return
elif [[ $# -lt 2 ]]; then
    echo "Merging all SAM files in this directory."
    samPrefix=""
else
    samPrefix=${2}
fi

finalFile=${1}

bamFiles=""

for file in $(find ${samPrefix}*.sam)
do
    echo "Converting file ${file}."
    samtools view -bS ${file} > ${file}.bam
    bamFiles="${bamFiles} ${file}.bam"
done

samtools merge ${finalFile} ${bamFiles}
