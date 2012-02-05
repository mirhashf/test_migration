#!/bin/bash -eu

# Fast Pipeline Runner
# Written by AJ Minich (aj.minich@binatechnologies.com), February 2012
#
# Runs the fast GATK pipeline on a sorted, marked BAM file.
# Logs output to <chrom>_<aligner>_<pipeline>.log
#

VERSION=1.01

# Constants
PICARD=~/programs/picard/dist
GATK=~/programs/gatk/dist
GATK_FAST=~/programs/gatk-fast
seqalto=/home/ajminich/programs/seq
bwa=/home/ajminich/programs/bwa/bwa-0.6.1/bwa
alignstats=/home/ajminich/programs/alignstats
SNP=/mnt/scratch0/public/data/variants/dbSNP/dbsnp_132.hg19.vcf

tmp=/mnt/scratch0/ajminich
pipeline_threads=1

if [[ $# -lt 4 ]]; then
    echo "Error: fewer than 4 params specified."
    exit
else
    reference=${1}
    reads=${2}
    aligner=${3}
    pipeline=${4}
    
    # Determine chromosome using the reference filename
    ref_file=$(basename $reference)
    chrom=${ref_file%.*}
    
    # Initialize the log file
    log_file="${chrom}_${aligner}_${pipeline}.log"
    
    echo "Full Pipeline Runner v${VERSION}" | tee ${log_file}
    echo "Aligner: ${aligner} with ${aligner_threads} threads" | tee -a ${log_file}
    echo "Pipeline: ${pipeline} with ${pipeline_threads} threads" | tee -a ${log_file}
    echo "Reference: ${reference}" | tee -a ${log_file}
    echo "Chromosome: ${chrom}" | tee -a ${log_file}
    echo "Reads: ${reads}" | tee -a ${log_file}
    echo "" | tee -a ${log_file}
    
fi

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

case "${aligner}" in

    "seqalto" | "seq" | "SeqAlto" | "Seqalto" )
        aligner="SeqAlto"
        ;;
        
    "bwa" | "BWA" )
        aligner="BWA"
        ;;
        
    * )
        echo "Error: I don't know how to use '${aligner}'."
        exit
        ;;

esac

aligned_reads=${reads}\_${aligner}\_AG\_sorted

echo -e "\n--------------------------- FAST GATK ---------------------------"

fast_gatk_start=$(timer)

java -Xms5g -Xmx5g -jar ${GATK_FAST}.jar \
    -T FastRealign \
    -I ${aligned_reads}\_marked.bam \
    -R ${reference} \
    -standard \
    -knownSites ${SNP} \
    -o {aligned_reads}\_realign.bam \
    -recalFile ${aligned_reads}\_realign.bam.csv
    #--suppress_recalibrate
    #--suppress_write_bam

echo "Fast GATK: $(timer ${fast_gatk_start})" | tee -a ${log_file}

echo -e "\n--------------------------- TABLE RECALIBRATION ---------------------------"

recalibrate_table_start=$(timer)

echo "Recalibrating base quality table."
if [ "`grep -v '#' ${aligned_reads}\_realign.bam.csv | grep -v "EOF" | wc -l`" = "1" ]
then
    cp ${aligned_reads}\_realign.bam ${aligned_reads}\_recalibrated.bam
else
    java -Xms5g -Xmx5g -jar ${GATK}/GenomeAnalysisTK.jar \
        -R ${reference} \
        -I ${aligned_reads}\_realign.bam  \
        -o ${aligned_reads}\_recalibrated.bam \
        -T TableRecalibration \
        -baq RECALCULATE \
        --doNotWriteOriginalQuals \
        -recalFile ${aligned_reads}\_realign.bam.csv \
        -et NO_ET
fi

# Rebuild index
java -Xms5g -Xmx5g -jar $PICARD/BuildBamIndex.jar INPUT=${aligned_reads}\_recalibrated.bam VALIDATION_STRINGENCY=LENIENT

echo "Table Recalibration: $(timer ${recalibrate_table_start})" | tee -a ${log_file}

echo -e "\n--------------------------- UNIFIED GENOTYPER VARIANT CALLING ---------------------------"

variant_calling_start=$(timer)

java -Xms5g -Xmx5g -Djava.io.tmpdir=${tmp} -jar $GATK/GenomeAnalysisTK.jar \
    -T UnifiedGenotyper \
    -I ${aligned_reads}\_recalibrated.bam \
    -R ${reference} \
    -D ${SNP} \
    -o ${aligned_reads}.vcf  \
    -et NO_ET \
    -dcov 1000 \
    -A AlleleBalance \
    -A DepthOfCoverage \
    -A MappingQualityZero \
    -baq CALCULATE_AS_NECESSARY \
    -stand_call_conf 30.0 \
    -stand_emit_conf 10.0 \
    -glm BOTH \
    -nt ${pipeline_threads}

echo "Variant Calling: $(timer ${variant_calling_start})" | tee -a ${log_file}

echo -e "\n--------------------------- PIPELINE PROCESSING COMPLETE ---------------------------"

echo "Total Running Time: $(timer ${start_time})"
echo "Log written to '${log_file}'."

