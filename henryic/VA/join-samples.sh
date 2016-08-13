#!/bin/bash

set -e

for SAMPLE in ${@}; do
    if [ ! -d ${SAMPLE} ]; then
        echo "Sample directory ${SAMPLE} not found"
        exit 1
    fi

    FLOWCELLS=$(find ${SAMPLE}/Fastq/ -mindepth 1 -maxdepth 1 -type d | cut -d"/" -f3)

    for FLOWCELL in ${FLOWCELLS}; do
        for LANE in $(ls ${SAMPLE}/Fastq/${FLOWCELL}/*.fastq.gz | cut -d "_" -f3 | sort | uniq); do
            for PAIR in 1 2;
                do cat ${SAMPLE}/Fastq/${FLOWCELL}/${SAMPLE}_NoIndex_${LANE}_R${PAIR}_*.fastq.gz > ${SAMPLE}/Fastq/${SAMPLE}_${FLOWCELL}_${LANE}_R${PAIR}.fastq.gz

                ./verify-cat.py ${SAMPLE}/Fastq/${SAMPLE}_${FLOWCELL}_${LANE}_R${PAIR}.fastq.gz ${SAMPLE}/Fastq/${FLOWCELL}/${SAMPLE}_NoIndex_${LANE}_R${PAIR}_*.fastq.gz
            done
        done
    done

    find ${SAMPLE}/Fastq/ -maxdepth 1 -type f -name "*.fastq.gz" -exec pigz -t {} \;

    shift

done

