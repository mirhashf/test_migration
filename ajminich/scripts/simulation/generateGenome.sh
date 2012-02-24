#! /bin/bash

# Genome Generator
# Written by AJ Minich (aj.minich@binatechnologies.com), January 2012
#
# This script generates a random genome. It should be run on
# a scratchdisk, since it outputs files to the local folder.

# Programs and Files
simuread=/home/ajminich/programs/loomis/software/loomis/simureads/simuRead
dbSNP=/mnt/scratch0/public/data/variants/dbSNP/CEU/CEU-1409-21.vcf
tmp=/mnt/scratch0/ajminich

# Local Variables
genome_sim=../simGenome

# Simulation Parameters
reference=/home/ajminich/data/hg19/ucsc.hg19.fasta     # run against full genome

pct_hetero=10               # genome heterozygosity
pct_random_mutation=0.15    # random mutation probability

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

echo -e "\n--------------------------- CREATING RANDOM GENOME ---------------------------"

simulation_start=$(timer)

# Create simulated genome
${simuread} random_genome \
    ${reference} \
    ${dbSNP} \
    ${genome_sim} \
    ${pct_hetero} \
    ${pct_random_mutation}

simulation_time=`echo $(timer ${simulation_start})`
echo "Random genome generation complete in ${simulation_time}."