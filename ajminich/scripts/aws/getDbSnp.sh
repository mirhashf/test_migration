#!/bin/bash

# Gets the hg19 dbSNP reference file.

DATA_FOLDER="/mnt/data/hg19/"

echo "Getting hg19 dbSNP reference data."
mkdir -p ${DATA_FOLDER}
s3cmd get --skip-existing s3://bina.data/hg19/dbsnp.hg19.vcf ${DATA_FOLDER}

echo "Completed downloading hg19 dbSNP reference data."
