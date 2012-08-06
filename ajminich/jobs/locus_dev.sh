#!/bin/bash

# Runs the Locus Dev datasets through Loomis.

SAMPLES="LS593 LS596 LS600 LS603 LS607"
FIRST_END="R1.fastq"
SECOND_END="R2.fastq"

BINA_DATA_PATH="bina://data/LocusDev/"
BINA_RESULTS_PATH="bina://jobs/LocusDev"

HOSTS="krakow-00:8080 krakow-01:8080"
KEY="gocardinal"

for sample in ${SAMPLES};
do
    echo "Queueing up Locus Dev sample ${sample}."
    echo bina --hosts ${HOSTS} --key ${KEY} run \
        -r ${BINA_DATA_PATH}/${sample}/${FIRST_END} \
        -R ${BINA_DATA_PATH}/${sample}/${SECOND_END} \
        -rg ${sample} -lb FMI -sm ${sample} \
        -o ${BINA_RESULTS_PATH} \
        -tg bina://data/LocusDev/SureSelect_038114_D_BED_20120106_merged_formatted.bed \
        --desc "Locus Dev Sample ${sample} (Local Mode)" -l
done

echo "All jobs queued."
