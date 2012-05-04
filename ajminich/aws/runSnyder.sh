#!/bin/bash

# Runs the Snyder EC2 script on the specified lanes, and uploads the results to EC2.
# Note: should be run as sudo.

SNYDER_EC2_SCRIPT=~/sandbox/ajminich/aws/seqaltoAlign.ec2.sh
S3_FOLDER=s3://bina.results/snyder/hg19.major/

bloodPrefix="A804NLABXX.s_"
bloodStart=6
bloodEnd=8

salivaPrefix="A808HKABXX.s_"
salivaStart=1
salivaEnd=8

function processAll() {
    
    prefix=${1}
    startIndex=${2}
    stopIndex=${3}
    
    for (( index=${startIndex}; index<=${stopIndex}; index++ ))
    do
    
        filename=${prefix}${index}
        echo "Processing Snyder lane ${filename}."

        ${SNYDER_EC2_SCRIPT} ${filename} 2>&1 | tee seqalto_${filename}.log

        # Upload the final BAM file and runlog to S3
        s3cmd put seqalto_${filename}.chr1.csorted.bam ${S3_FOLDER}
        s3cmd put seqalto_${filename}.log ${S3_FOLDER}
    
        echo "Snyder lane ${filename} complete."
        
    done

}

processAll ${bloodPrefix} ${bloodStart} ${bloodEnd}
processAll ${salivaPrefix} ${salivaStart} ${salivaEnd}
