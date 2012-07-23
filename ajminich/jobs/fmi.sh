#!/bin/bash

# Runs the FMI datasets through Loomis.

FMI_DATA_FOLDER="/mnt/scratch0/public/data/FMI"
FIRST_END_SUFFIX=".read1.fq"
SECOND_END_SUFFIX=".read2.fq"

# Get list of prefixes to process
files=`ls ${FMI_DATA_FOLDER} | cut -f1,2 -d'.' | uniq`

echo ${files}