#!/bin/bash

# Runs the FMI datasets through Loomis.

FMI_DATA_FOLDER="/export/data/FMI/P0085/reads"
FIRST_END_SUFFIX="read1.fastq"
SECOND_END_SUFFIX="read2.fastq"

BINA_DATA_PATH="bina://data/FMI/P0085/reads"
BINA_RESULTS_PATH="bina://jobs/FMI"

HOSTS="krakow-00:8080 krakow-01:8080"
KEY="gocardinal"

# Get list of prefixes to process
files=`bina ls bina:///data/FMI/P0085/reads | grep \.fastq | cut -d'/' -f8 | cut -d'.' -f1-3 | uniq`
files=`echo $files | fmt -1 | head -1`

for file in ${files};
do
    sample_name=`basename ${file}`

    echo "Queueing up FMI sample ${sample_name}."
    bina --hosts ${HOSTS} --key ${KEY} run \
        -r ${BINA_DATA_PATH}/${sample_name}.${FIRST_END_SUFFIX} \
        -R ${BINA_DATA_PATH}/${sample_name}.${SECOND_END_SUFFIX} \
        -rg ${sample_name} -lb FMI -sm ${sample_name} \
        -o ${BINA_RESULTS_PATH} --use_bwa \
        --no_genotyper \
        --debug \
        --address napa@binatechnologies.com --sender_id binabox \
        --desc "FMI Sample ${sample_name} (Local Mode)" -l
done

echo "All jobs queued."

