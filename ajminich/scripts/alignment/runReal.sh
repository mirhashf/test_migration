#! /bin/bash

# Real Runner
# Written by AJ Minich (aj.minich@binatechnologies.com), January 2012
#
# This script runs a single run-through on the given real data. It should be run
# on a scratchdisk, since it outputs files to the local folder.

# Programs and Files
seqalto=/home/ajminich/programs/seq
bwa=/home/ajminich/programs/bwa/bwa-0.6.1/bwa
trimmer=/home/ajminich/programs/simulation/trimBWAstyle.pl
alignstats=/home/ajminich/programs/alignstats
tmp=/mnt/scratch0/ajminich

# Local Variables
read_names=A804NLABXX\_lane1
seqalto_align=seqAlign
bwa_align=bwaAlign
results_folder=results
results_file=real_reads

# Parameters
reference=/home/ajminich/data/hg19/ucsc.hg19.fasta
reads1=/mnt/scratch0/public/data/snyder/A804NLABXX.s\_1\_1.fq.gz
reads2=/mnt/scratch0/public/data/snyder/A804NLABXX.s\_1\_2.fq.gz

num_pairs=100000          # total number of read pairs
quality=20                  # Mapping quality threshold for read trimming

# Aligner Options
num_threads=16

# correctness
wiggle=21
num_pairs_half=`expr ${num_pairs} / 2`
min_mapq=1
min_len=1

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

#echo -e "\n--------------------------- CREATING REAL READS ---------------------------"

#subsampling_start=$(timer)

# Subsample reads
#echo "Subsampling ${reads1} and ${reads2} down to ${num_pairs} read pairs."
#${seqalto} sample_reads \
#    ${num_pairs} \
#    ${read_names} \
#    ${reads1} \
#    ${reads2}

# Trim reads
#cat ${read_names}\_1.fq | perl ${trimmer} -q ${quality} > ${read_names}\_trimmed\_1.fq
#cat ${read_names}\_2.fq | perl ${trimmer} -q ${quality} > ${read_names}\_trimmed\_2.fq

#subsampling_time=`echo $(timer ${subsampling_start})`
#echo "Read subsampling complete in ${subsampling_time}."

echo -e "\n--------------------------- RUNNING SEQALTO ---------------------------"

seqalto_start=$(timer)

# Run SeqAlto
echo "Running SeqAlto aligner."
time ${seqalto} align \
    --idx ${reference}\_22.sidx \
    -p ${num_threads} \
    -1 ${read_names}\_trimmed\_1.fq \
    -2 ${read_names}\_trimmed\_2.fq \
    > ${seqalto_align}.sam

seqalto_time=`echo $(timer ${seqalto_start})`
echo "SeqAlto alignment complete in ${seqalto_time}."

# Convert to BAM
echo "Converting SeqAlto SAM results file to BAM."
samtools view -bS ${seqalto_align}.sam > ${seqalto_align}.bam

echo -e "\n--------------------------- RUNNING BWA ---------------------------"

bwa_start=$(timer)

# Run BWA
echo "Running BWA aligner."
time ${bwa} aln \
    ${reference} \
    ${read_names}\_trimmed\_1.fq \
    -t ${num_threads} \
    -f ${read_names}\_1.sai
time ${bwa} aln \
    ${reference} \
    ${read_names}\_trimmed\_2.fq \
    -t ${num_threads} \
    -f ${read_names}\_2.sai
    #-q ${quality} \

echo "Performing BWA paired-end merging."
time ${bwa} sampe \
    -P ${reference} \
    ${read_names}\_1.sai \
    ${read_names}\_2.sai \
    ${read_names}\_trimmed\_1.fq \
    ${read_names}\_trimmed\_2.fq \
    > ${bwa_align}.sam

bwa_time=`echo $(timer ${bwa_start})`
echo "BWA alignment complete in ${bwa_time}."

# Convert to BAM
echo "Converting BWA SAM results file to BAM."
samtools view -bS ${bwa_align}.sam > ${bwa_align}.bam

echo -e "\n--------------------------- RUNNING ALIGNSTATS ---------------------------"

alignstats_start=$(timer)

# Run AlignDiff
java -jar ${alignstats} \
    ${seqalto_align}.bam ${bwa_align}.bam \
    ${results_folder}/ \
    --pivot-report ${results_file} \
    --temp-dir ${tmp}
#    --diff \
#    --output-html ${results_file}.html \

alignstats_time=`echo $(timer ${alignstats_start})`
echo "AlignStats analysis complete in ${alignstats_time}."

echo "Alignment Difference report written to '${results_file}'."

echo -e "\n--------------------------- SIMULATION COMPLETE ---------------------------\n"

total_time=`echo $(timer ${start_time})`

echo -e "Simulation Generation: \t ${simulation_time}"
echo -e "SeqAlto Alignment: \t ${seqalto_time}"
echo -e "BWA Alignment: \t ${bwa_time}"
echo -e "AlignStats Analysis: \t ${alignstats_time}"
echo -e "Total Running Time: \t ${total_time}"