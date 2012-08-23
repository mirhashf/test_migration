#!/bin/bash

# Runs the FMI datasets through Loomis.

FMI_DATA_FOLDER="/export/krakow/data/FMI/P0085/reads"
FIRST_END_SUFFIX="read1.fastq"
SECOND_END_SUFFIX="read2.fastq"

BINA_DATA_PATH="bina://data/FMI/P0085/reads"
BINA_RESULTS_PATH="bina://jobs/FMI"

HOSTS="krakow-00:8080 krakow-01:8080"
KEY="gocardinal"

# Get list of prefixes to process
files=`ls ${FMI_DATA_FOLDER}/*.fastq | cut -f1-3 -d'.' | uniq`

for file in ${files};
do
    sample_name=`basename ${file}`
    
    echo "Queueing up FMI sample ${sample_name}."
    /home/ajminich/sdk/bin/bina --hosts ${HOSTS} --key ${KEY} run \
        -r ${BINA_DATA_PATH}/${sample_name}.${FIRST_END_SUFFIX} \
        -R ${BINA_DATA_PATH}/${sample_name}.${SECOND_END_SUFFIX} \
        -rg ${sample_name} -lb FMI -sm ${sample_name} \
        -o ${BINA_RESULTS_PATH} --use_bwa \
        --desc "FMI Sample ${sample_name} (Local Mode)" -l
done

echo "All jobs queued."
