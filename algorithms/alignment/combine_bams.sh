#!/bin/bash -eu

port1=1367
port2=1368
machine_name="tehran"
machine_nums="00 01 02 03 06"
num_threads=24
gatk="/home/ajminich/seqalto/third-party/gatk/dist/GenomeAnalysisTK.jar"
ref="/home/ajminich/genome/human/CEUref/CEUref.fasta"
bed_file="/export/tehran/data/snyder/pcr/cell_6126_mmc1.bed"

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
    for thread_num in $(seq 0 `expr ${num_threads} - 1`)
    do
        bam_files="${bam_files} -I ${folder}/${thread_num}.bam"
    done
done

echo "Initiating merge."
java -Xms15g -Xmx15g -jar ${gatk} -T PrintReads \
    ${bam_files} \
    --out ${final_file} \
    -R ${ref} \
    --validation_strictness LENIENT \
    -L ${bed_file}