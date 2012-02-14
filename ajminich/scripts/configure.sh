#! /bin/bash

# populateData.sh
# Populates data on the scratchdisk using s3cmd.

EXEC_FOLDER="/mnt/execution"

echo "Acquiring hg19 genome data from S3."

INDEX_LOC="s3://bina.data/hg19"
INDEX_FILE="ucsc.hg19.fasta"
INDEX_FOLDER="/mnt/data/hg19"

mkdir ${INDEX_FOLDER}

s3cmd get ${INDEX_LOC}/${INDEX_FILE} ${INDEX_FOLDER}/${INDEX_FILE}
s3cmd get ${INDEX_LOC}/${INDEX_FILE}.ann ${INDEX_FOLDER}/${INDEX_FILE}.ann
s3cmd get ${INDEX_LOC}/${INDEX_FILE}.amb ${INDEX_FOLDER}/${INDEX_FILE}.amb
s3cmd get ${INDEX_LOC}/${INDEX_FILE}.bwt ${INDEX_FOLDER}/${INDEX_FILE}.bwt
s3cmd get ${INDEX_LOC}/${INDEX_FILE}.rbwt ${INDEX_FOLDER}/${INDEX_FILE}.rbwt
s3cmd get ${INDEX_LOC}/${INDEX_FILE}.pac ${INDEX_FOLDER}/${INDEX_FILE}.pac
s3cmd get ${INDEX_LOC}/${INDEX_FILE}.rpac ${INDEX_FOLDER}/${INDEX_FILE}.rpac
s3cmd get ${INDEX_LOC}/${INDEX_FILE}.sa ${INDEX_FOLDER}/${INDEX_FILE}.sa
s3cmd get ${INDEX_LOC}/${INDEX_FILE}.rsa ${INDEX_FOLDER}/${INDEX_FILE}.rsa

echo "Acquiring reads data from S3."

READS_LOC="s3://bina.data/reads"
READS_FILE="sim_hg19_20M"
READS_FOLDER="/mnt/data/reads"

mkdir ${READS_FOLDER}

# Get the data
s3cmd get ${READS_LOC}/${READS_FILE}\_1.fq ${READS_FOLDER}/${READS_FILE}\_1.fq
s3cmd get ${READS_LOC}/${READS_FILE}\_2.fq ${READS_FOLDER}/${READS_FILE}\_2.fq

echo "Finished populating data."

mkdir ${EXEC_FOLDER}
chmod a+rw ${EXEC_FOLDER}
