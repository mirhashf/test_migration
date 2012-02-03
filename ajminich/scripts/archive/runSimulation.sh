#! /bin/bash

# Simulation Runner
# Written by AJ Minich (aj.minich@binatechnologies.com), January 2012
#
# This script runs a single simulation on the given data. It should be run on
# a scratchdisk, since it outputs files to the local folder.

# Programs and Files
simuread=/home/ajminich/programs/loomis/software/loomis/simureads/simuRead
dbSNP=/mnt/scratch0/public/data/variants/dbSNP/CEU/CEU-1409-21.vcf
seqalto=/home/ajminich/programs/seq
bwa=/home/ajminich/programs/bwa/bwa-0.6.1/bwa
trimmer=/home/ajminich/programs/align/trimBWAstyle.pl
alignstats=/home/ajminich/programs/alignstats
tmp=/mnt/scratch0/ajminich

# Local Variables
genome_sim=../simGenome
read_names=reads
seqalto_align=seqAlign
bwa_align=bwaAlign

# Simulation Parameters
reference=/home/ajminich/data/hg19/ucsc.hg19.fasta     # run against full genome

pct_hetero=10               # genome heterozygosity
pct_random_mutation=0.15    # random mutation probability
read_length_1=100           # length of 5' reads
read_length_2=100           # length of 3' reads
num_pairs=100000            # total number of read pairs
length_avg=330              # average of template lengths
length_std=7                # standard deviation of template lengths
pct_unmappable=3            # percentage of truly random read data

# Output parameters
results_folder=results
results_file=sim_reads_${num_pairs}

# Aligner Options
num_threads=16
max_template_size=400       # Max template size (SeqAlto)
quality=20                  # Mapping quality threshold for read trimming

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

echo -e "\n--------------------------- CREATING SIMULATED READS ---------------------------"

simulation_start=$(timer)

# Create simulated genome
#echo "Creating simulated genome."
#${simuread} random_genome \
#    ${reference} \
#    ${dbSNP} \
#    ${genome_sim} \
#    ${pct_hetero} \
#    ${pct_random_mutation}

# Create simulated reads
#echo "Creating simulated reads."
#${simuread} pair_reads \
#    ${genome_sim} \
#    ${read_length_1} \
#    ${read_length_2} \
#    ${num_pairs} \
#    ${length_avg} \
#    ${length_std} \
#    ${pct_unmappable} \
#    ${read_names}

# Trim reads
#cat ${read_names}\_1.fq | perl ${trimmer} -q ${quality} > ${read_names}\_trimmed\_1.fq
#cat ${read_names}\_2.fq | perl ${trimmer} -q ${quality} > ${read_names}\_trimmed\_2.fq
#
#simulation_time=`echo $(timer ${simulation_start})`
#echo "Simulated read generation complete in ${simulation_time}."

echo -e "\n--------------------------- RUNNING SEQALTO ---------------------------"

seqalto_load_start=$(timer)

echo "Loading the genome reference '${reference}'."
#nohup ${seqalto} load_index ${reference}\_22.sidx &

seqalto_load_time=`echo $(timer ${seqalto_load_start})`
echo "SeqAlto reference genome index loading complete in ${seqalto_load_time}."

seqalto_align_start=$(timer)

# Run SeqAlto
echo "Running SeqAlto aligner."
time ${seqalto} align \
    --idx ${reference}\_22.sidx \
    -m ${length_avg} \
    -i ${max_template_size} \
    -p ${num_threads} \
    -1 ${read_names}\_trimmed\_1.fq \
    -2 ${read_names}\_trimmed\_2.fq \
    > ${seqalto_align}.sam
    #--trim ${quality}
    
seqalto_align_time=`echo $(timer ${seqalto_align_start})`
echo "SeqAlto alignment complete in ${seqalto_align_time}."

# Perform pair-alignment checking
echo "Checking pairs."
grep "pa;" ${seqalto_align}.sam > ${seqalto_align}.pa.sam
grep "ma;" ${seqalto_align}.sam > ${seqalto_align}.ma.sam

${seqalto} check_pair \
    ${genome_sim}.pa.map \
    ${seqalto_align}.pa.sam \
    ${wiggle} \
    ${num_pairs_half} \
    ${min_mapq} \
    ${min_len} \
    ${seqalto_align}\_checkPa \
    > ${seqalto_align}\_pa.out
${seqalto} check_pair \
    ${genome_sim}.ma.map \
    ${seqalto_align}.ma.sam \
    ${wiggle} \
    ${num_pairs_half} \
    ${min_mapq} \
    ${min_len} \
    ${seqalto_align}\_checkMa \
    > ${seqalto_align}\_ma.out

# Convert to BAM
echo "Converting SAM results file to BAM."
samtools view -bS ${seqalto_align}.sam > ${seqalto_align}.bam

echo -e "\n--------------------------- RUNNING BWA ---------------------------"

bwa_start=$(timer)

# Run BWA
echo "Running BWA aligner."
time ${bwa} aln \
    ${reference} \
    ${read_names}\_trimmed\_1.fq \
    -t ${num_threads} \
    > ${read_names}\_1.sai
time ${bwa} aln \
    ${reference} \
    ${read_names}\_trimmed\_2.fq \
    -t ${num_threads} \
    > ${read_names}\_2.sai
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

# Perform pair-alignment checking
echo "Checking pairs."
grep "pa;" ${bwa_align}.sam > ${bwa_align}.pa.sam
grep "ma;" ${bwa_align}.sam > ${bwa_align}.ma.sam

${seqalto} check_pair \
    ${genome_sim}.pa.map \
    ${bwa_align}.pa.sam \
    ${wiggle} \
    ${num_pairs_half} \
    ${min_mapq} \
    ${min_len} \
    ${bwa_align}\_checkPa \
    > ${bwa_align}\_pa.out
${seqalto} check_pair \
    ${genome_sim}.ma.map \
    ${bwa_align}.ma.sam \
    ${wiggle} \
    ${num_pairs_half} \
    ${min_mapq} \
    ${min_len} \
    ${bwa_align}\_checkMa \
    > ${bwa_align}\_ma.out

# Convert to BAM
echo "Converting SAM results file to BAM."
samtools view -bS ${bwa_align}.sam > ${bwa_align}.bam

echo -e "\n--------------------------- RUNNING ALIGNSTATS ---------------------------"

alignstats_start=$(timer)

# Run AlignDiff
java -jar ${alignstats} \
    ${seqalto_align}.bam ${bwa_align}.bam \
    ${results_folder}/ \
    --pivot-report ${results_file} \
    --report ${results_file}.html \
    --temp-dir ${tmp}
    
alignstats_time=`echo $(timer ${alignstats_start})`
echo "AlignStats analysis complete in ${alignstats_time}."

echo "Alignment Difference report written to '${results_file}'."

echo -e "\n--------------------------- SIMULATION COMPLETE ---------------------------"

total_time=`echo $(timer ${start_time})`

echo -e "Simulation Generation: \t ${simulation_time}"
echo -e "SeqAlto Index Loading: \t ${seqalto_load_time}"
echo -e "SeqAlto Alignment: \t ${seqalto_align_time}"
echo -e "BWA Alignment: \t ${bwa_time}"
echo -e "AlignStats Analysis: \t ${alignstats_time}"
echo -e "Total Running Time: \t ${total_time}"