#! /bin/bash

# Reads Generator
# Written by AJ Minich (aj.minich@binatechnologies.com), January 2012
#
# This script generates simulated reads for the given genome. It should be run on
# a scratchdisk, since it outputs files to the local folder.

num_pairs=${1}

# Programs and Files
simuread=/home/ajminich/programs/simuread
trimmer=/home/ajminich/programs/align/trimBWAstyle.pl
tmp=/mnt/scratch0/ajminich

# Local Variables
genome_sim=../simGenome
read_names=reads

# Simulation Parameters
reference=/home/ajminich/data/hg19/ucsc.hg19.fasta     # run against full genome

read_length_1=100           # length of 5' reads
read_length_2=100           # length of 3' reads
length_avg=330              # average of template lengths
length_std=7                # standard deviation of template lengths
pct_unmappable=3            # percentage of truly random read data

quality=20

# Timer for determining elapsed time of each operation.
function timer()
{
    if [[ $# -eq 0 ]]; then
        echo $(date '+%s')
    else
        local  stime=$1
        etime=$(date '+%s')

        if [[ -z "$stime" ]]; then stime=$etime; fi

        dt=$((etime - stime))
        ds=$((dt % 60))
        dm=$(((dt / 60) % 60))
        dh=$((dt / 3600))
        printf '%d:%02d:%02d' $dh $dm $ds
    fi
}

echo -e "\n--------------------------- CREATING SIMULATED READS ---------------------------"

simulation_start=$(timer)

# Create simulated reads
echo "Creating simulated reads."
${simuread} pair_reads \
    ${genome_sim} \
    ${read_length_1} \
    ${read_length_2} \
    ${num_pairs} \
    ${length_avg} \
    ${length_std} \
    ${pct_unmappable} \
    ${read_names}

# Trim reads
cat ${read_names}\_1.fq | perl ${trimmer} -q ${quality} > ${read_names}\_trimmed\_1.fq
cat ${read_names}\_2.fq | perl ${trimmer} -q ${quality} > ${read_names}\_trimmed\_2.fq

simulation_time=`echo $(timer ${simulation_start})`
echo "Simulated read generation complete in ${simulation_time}."
