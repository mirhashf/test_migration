#!/bin/bash

# Looks through each node and determines which node has your job

if [ -z ${1} ] || [ -z ${2} ]; then
    echo "Usage: ${0} <box_name> <job_id_prefix>"
    exit
fi

box=${1}
id=${2}

for i in `seq -s " " -f %02g 0 3` 
do
     echo "${box}-${i}:"
     ssh binatech@${box}-${i} ls /mnt/scratch/local/jobs/ | grep "${id}*"
done
