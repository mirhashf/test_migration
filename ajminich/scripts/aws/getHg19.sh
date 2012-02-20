#!/bin/bash

# Gets the reference data for all of hg19.

DATA_FOLDER="/mnt/data/hg19/"

echo "Getting hg19 reference data."
mkdir -p ${DATA_FOLDER}
s3cmd get s3://bina.data/hg19/ucsc.hg19.* ${DATA_FOLDER}

echo "Completed downloading hg19 reference data."