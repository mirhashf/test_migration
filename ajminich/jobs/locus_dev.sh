#!/bin/bash

# Runs the Locus Dev datasets through Loomis.

SAMPLES="LS593 LS596 LS600 LS603 LS607"
FIRST_END="R1.fastq"
SECOND_END="R2.fastq"

BINA_DATA_PATH="bina://data/locus_dev/data"
BINA_RESULTS_PATH="bina://jobs/LocusDev"

HOSTS="milan-00:8080 milan-01:8080"
KEY="gocardinal"

for sample in ${SAMPLES};
do
    echo "Queueing up Locus Dev sample ${sample}."
    echo bina --hosts ${HOSTS} --key ${KEY} run \
        -r ${BINA_DATA_PATH}/${sample}/${FIRST_END} \
        -R ${BINA_DATA_PATH}/${sample}/${SECOND_END} \
        -rg ${sample} -lb FMI -sm ${sample} \
        -o ${BINA_RESULTS_PATH} \
        --desc "Locus Dev Sample ${sample}" \
        --enable_local_run
done

echo "All jobs queued."
