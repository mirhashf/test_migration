#!/bin/bash -eu

# Remove FASTQ Additional Info
# Written by AJ Minich (aj.minich@binatechnologies.com) with
# Josh Paul (jpaul@binatechnologies.com), March 2012
#
# Strips the additional reads info from FASTQ files.
# Essentially, any line that begins with + will be replaced
# with a line only containing '+'.

if [[ $# -lt 1 ]]; then
    echo "FASTQ Additional Info Removal"
    echo "Usage: remove_fastq_info.sh <read_prefix>"
    echo ""
    echo "Removes additional reads info on lines beginning with '+'."
else
    sed 's/^+.*/+/' ${1}_1.fq > ${1}_cleaned_1.fq
    sed 's/^+.*/+/' ${1}_2.fq > ${1}_cleaned_2.fq
fi