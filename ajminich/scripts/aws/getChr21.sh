#!/bin/bash

# Gets the hg19 reference data for Chromosome 21.

DATA_FOLDER="/mnt/data/chr21/"

echo "Getting hg19 chr21 reference data."
mkdir -p ${DATA_FOLDER}
s3cmd get s3://bina.data/hg19/chr21.fa* ${DATA_FOLDER}

echo "Completed downloading chr21 reference data."