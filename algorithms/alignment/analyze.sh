#!/bin/bash -eu

# Analyzer
# Written by AJ Minich (aj.minich@binatechnologies.com), February 2012
#
# Performs check_pair analysis on the provided file.

seqalto="/home/ajminich/programs/seq"

if [[ $# -lt 3 ]]; then
    echo "Simulated Data Alignment Analyzer"
    echo "Usage: analyze.sh <genome_prefix> <num_pairs> <alignment_file>.sam"
    echo ""
    echo "Performs check_pairs analysis on the provided alignment."
else
    genome_sim=${1}
    num_pairs=${2}
    align=${3}
    
    wiggle=21
    num_pairs_half=`expr ${num_pairs} / 2`
    min_mapq=1
    min_len=1
    
    echo "Checking pairs."
    grep "pa;" ${align}.sam > ${align}.pa.sam
    grep "ma;" ${align}.sam > ${align}.ma.sam

    ${seqalto} -mode check_pair \
        -map_file ${genome_sim}.pa.map \
        -sam_file ${align}.pa.sam \
        -wiggle_bp ${wiggle} \
        -num_pairs ${num_pairs_half} \
        -min_mapq ${min_mapq} \
        -min_len ${min_len} \
        -output_prefix ${align}\_checkPa \
        > ${align}\_pa.out
    ${seqalto} -mode check_pair \
        -map_file ${genome_sim}.ma.map \
        -sam_file ${align}.ma.sam \
        -wiggle_bp ${wiggle} \
        -num_pairs ${num_pairs_half} \
        -min_mapq ${min_mapq} \
        -min_len ${min_len} \
        -output_prefix ${align}\_checkMa \
        > ${align}\_ma.out
fi