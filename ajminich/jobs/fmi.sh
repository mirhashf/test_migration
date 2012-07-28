#!/bin/bash

# Runs the FMI datasets through Loomis.

FMI_DATA_FOLDER="/mnt/scratch0/public/data/FMI"
FIRST_END_SUFFIX="read1.fq"
SECOND_END_SUFFIX="read2.fq"

BINA_DATA_PATH="bina://data/FMI"
BINA_RESULTS_PATH="bina://jobs/FMI"

HOSTS="milan-00:8080 milan-01:8080"
KEY="gocardinal"

# Get list of prefixes to process
files=`ls ${FMI_DATA_FOLDER}/*.fq | cut -f1,2 -d'.' | uniq`

for file in ${files};
do
    sample_name=`basename ${file}`
    
    echo "Queueing up FMI sample ${sample_name}."
    bina --hosts ${HOSTS} --key ${KEY} run \
        -r ${BINA_DATA_PATH}/${sample_name}.${FIRST_END_SUFFIX} \
        -R ${BINA_DATA_PATH}/${sample_name}.${SECOND_END_SUFFIX} \
        -rg ${sample_name} -lb FMI -sm ${sample_name} \
        -o ${BINA_RESULTS_PATH} --use_bwa --use_broad_gatk \
        --desc "FMI Sample ${sample_name} using Broad GATK"
done

echo "All jobs queued."
