#! /bin/bash

# Comparison Runner
# Written by AJ Minich (aj.minich@binatechnologies.com), January 2012
#
# This script runs a single comparison of the given data. It should be run on
# a scratchdisk, since it outputs files to the local folder.

# Inputs
read_names=${1}
prefix=${2}

# Programs and Files
seqalto=/home/ajminich/programs/seq
bwa=/home/ajminich/programs/bwa
alignstats=/home/ajminich/programs/alignstats
tmp=/mnt/scratch0/ajminich

# Local Variables
seqalto_align=seqAlign
bwa_align=bwaAlign

# Simulation Parameters
reference=/home/ajminich/data/hg19/ucsc.hg19.fasta     # run against full genome
#reference=/home/ajminich/data/hg19/chr21.fa

# Output parameters
results_folder=${prefix}
results_file=${prefix}_results

# Aligner Options
num_threads=16
quality=20                  # Mapping quality threshold for read trimming

# Timer for determining elapsed time of each operation.
function timer()
{
    if [[ $# -eq 0 ]]; then
        echo $(date '+%s')
    else
        local  stime=$1
        etime=$(date '+%s')

        if [[ -z "$stime" ]]; then stime=$etime; fi

        dt=$((etime - stime))
        ds=$((dt % 60))
        dm=$(((dt / 60) % 60))
        dh=$((dt / 3600))
        printf '%d:%02d:%02d' $dh $dm $ds
    fi
}

start_time=$(timer)

echo -e "\n--------------------------- RUNNING SEQALTO ---------------------------"

seqalto_align_start=$(timer)

# Run SeqAlto
echo "Running SeqAlto aligner."
time ${seqalto} align \
    --idx ${reference}\_22.midx \
    -p ${num_threads} \
    -1 ${read_names}\_1.fq \
    -2 ${read_names}\_2.fq \
    > ${seqalto_align}.sam
    #--trim ${quality}
    
seqalto_align_time=`echo $(timer ${seqalto_align_start})`
echo "SeqAlto alignment complete in ${seqalto_align_time}."

# Convert to BAM
seqalto_bam_conversion_start=$(timer)

echo "Converting SAM results file to BAM."
samtools view -bS ${seqalto_align}.sam > ${seqalto_align}.bam

seqalto_bam_conversion_time=`echo $(timer ${seqalto_bam_conversion_start})`
echo "SeqAlto SAM-to-BAM conversion complete in ${seqalto_bam_conversion_time}."

echo -e "\n--------------------------- RUNNING BWA ---------------------------"

bwa_start=$(timer)

# Run BWA
echo "Running BWA aligner."
time ${bwa} aln \
    ${reference} \
    ${read_names}\_1.fq \
    -t ${num_threads} \
    > ${read_names}\_1.sai
time ${bwa} aln \
    ${reference} \
    ${read_names}\_2.fq \
    -t ${num_threads} \
    > ${read_names}\_2.sai
    #-q ${quality} \

echo "Performing BWA paired-end merging."
time ${bwa} sampe \
    -P ${reference} \
    ${read_names}\_1.sai \
    ${read_names}\_2.sai \
    ${read_names}\_1.fq \
    ${read_names}\_2.fq \
    > ${bwa_align}.sam

bwa_time=`echo $(timer ${bwa_start})`
echo "BWA alignment complete in ${bwa_time}."

# Convert to BAM
bwa_bam_conversion_start=$(timer)

echo "Converting SAM results file to BAM."
samtools view -bS ${bwa_align}.sam > ${bwa_align}.bam

bwa_bam_conversion_time=`echo $(timer ${bwa_bam_conversion_start})`
echo "BWA SAM-to-BAM conversion complete in ${bwa_bam_conversion_time}."

#echo -e "\n--------------------------- RUNNING ALIGNSTATS ---------------------------"
#
#alignstats_start=$(timer)
#
## Run AlignDiff
#java -jar ${alignstats} \
#    ${seqalto_align}.bam ${bwa_align}.bam \
#    ${results_folder}/ \
#    --report ${results_file}.html \
#    --temp-dir ${tmp}
##    --pivot-report ${results_file} \
#    
#alignstats_time=`echo $(timer ${alignstats_start})`
#echo "AlignStats analysis complete in ${alignstats_time}."
#
#echo "Alignment Difference report written to '${results_file}'."

echo -e "\n--------------------------- SIMULATION COMPLETE ---------------------------"

total_time=`echo $(timer ${start_time})`

echo -e "SeqAlto Alignment:      \t ${seqalto_align_time}"
echo -e "SAM-to-Bam Conversion:  \t ${seqalto_bam_conversion_time}"
echo -e "BWA Alignment:          \t ${bwa_time}"
echo -e "SAM-to-Bam Conversion:  \t ${bwa_bam_conversion_time}"
#echo -e "AlignStats Analysis: \t ${alignstats_time}"
echo -e "Total Running Time:     \t ${total_time}"