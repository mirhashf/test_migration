#!/bin/bash -eu

port1=1367
port2=1368
machine_name="tehran"
machine_nums="00 01 02 03 06"
picard_folder="/home/ajminich/programs/picard/dist"

if [[ $# -lt 1 ]]; then
    echo "Bina Box BAM Combiner"
    echo "Merges BAM files in the Bina Box output folder into a single BAM file."
    echo "Usage: combine_bams.sh [final_file.bam]"
    echo "  <final_file.bam>   - the final BAM file to output"
    exit
fi

final_file=${1}

bam_files=""

for machine_num in ${machine_nums}
do
    folder="${machine_name}-${machine_num}-${port1}-${port2}"
    for file in $(find ${folder}/*.bam)
    do
        echo "Discovered file ${file}."
        bam_files="${bam_files} I=${file}"
    done
done

echo "Initiating merge."
java -Xms15g -Xmx15g -jar ${picard_folder}/MergeSamFiles.jar ${bam_files} \
    OUTPUT=${final_file} \
    SORT_ORDER=coordinate \
    VALIDATION_STRINGENCY=LENIENT \
    USE_THREADING=T MSD=true