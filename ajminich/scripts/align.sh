#! /bin/bash

# Alignment Runner
# Written by AJ Minich (aj.minich@binatechnologies.com), February 2012
#
# This script runs an alignment of the given data using the selected aligners.
# It should be run on a scratchdisk, since it outputs files to the local folder.

VERSION=1.00

# Programs and Files
seqalto=/home/ajminich/programs/seq
bwa=/home/ajminich/programs/bwa
alignstats=/home/ajminich/programs/alignstats
tmp=/mnt/scratch0/ajminich

# Aligner Options
kmer_size=22
num_threads=16
quality=20                  # Mapping quality threshold for read trimming

if [[ $# -lt 4 ]]; then
    echo "Alignment Runner"
    echo "Runs alignment on the specified reads files."
    echo ""
    echo "Usage: align.sh <reference> <reads_prefix> <aligner> <run_alignstats>"
    echo ""
    echo "Reference:            location of the reference FASTA file, pre-indexed"
    echo "Reads prefix:         reads will be taken from <reads_prefix>_1.fq and <reads_prefix>_2.fq"
    echo "Aligner choices:      seqalto, bwa, both, none"
    echo "Alignstats options:   1/yes/Y/y to run, anything else to not run."
    echo ""
    echo "Notes:"
    echo "- Alignment is run with ${num_threads} thread(s)."
    echo "- Assuming k-mer size of ${kmer_size} for SeqAlto index."
    #echo "- Reads are trimmed with a quality threshold of ${quality}."
    exit
else
    reference=${1}
    reads_prefix=${2}
    aligner=${3}
    runstats=${4}
    
    # Initialize the log file
    log_file="${reads_prefix}_${aligner}.log"
    
    echo "Alignment Runner v${VERSION}" | tee ${log_file}
    echo "Aligner: ${aligner} with ${num_threads} thread(s)" | tee -a ${log_file}
    echo "Reference: ${reference}" | tee -a ${log_file}
    echo "Reads: ${reads_prefix}" | tee -a ${log_file}
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

    "seqalto" | "seq" | "SeqAlto" | "Seqalto" | "both" | "all" | "BOTH" | "ALL" )
        
        seqalto_start=$(timer)
    
        echo -e "\n--------------------------- ALIGNMENT WITH SEQALTO ---------------------------"

        #index_load_start=$(timer)
        #${seqalto} load_index ${reference}\_22.midx
        #echo "SeqAlto Index Loading: $(timer ${index_load_start})" | tee -a ${log_file}

        ${seqalto} align \
            --idx ${reference}_${kmer_size}.midx \
            -p ${num_threads} \
            -1 ${reads_prefix}_1.fq \
            -2 ${reads_prefix}_2.fq \
            > ${reads}\_SeqAlto.sam
            
        echo "SeqAlto Alignment: $(timer ${seqalto_start})" | tee -a ${log_file}
            
        ;;
        
    "bwa" | "BWA" | "both" | "all" | "BOTH" | "ALL"  )
    
        bwa_start=$(timer)
    
        echo -e "\n--------------------------- ALIGNMENT WITH BWA ---------------------------"

        ${bwa} aln \
            ${reference} \
            ${reads_prefix}_1.fq \
            -t ${num_threads} \
            > ${reads_prefix}_1.sai
        ${bwa} aln \
            ${reference} \
            ${reads_prefix}_2.fq \
            -t ${num_threads} \
            > ${reads_prefix}_2.sai
        
        echo "Performing BWA paired-end merging."
        ${bwa} sampe \
            -P ${reference} \
            ${reads_prefix}_1.sai \
            ${reads_prefix}_2.sai \
            ${reads_prefix}_1.fq \
            ${reads_prefix}_2.fq \
            > ${reads_prefix}_BWA.sam
        
        echo "BWA Alignment: $(timer ${bwa_start})" | tee -a ${log_file}
        
        ;;
        
    "none" | "neither" )
    
        echo "No alignments performed."
        
        ;;
        
    * )
        echo "Error: I don't know how to use '${aligner}'." | tee -a ${log_file}
        exit
        ;;

esac

# AlignStats
case "${runstats}" in

    "1" | "yes" | "Y" | "y" )
    
        echo -e "\n--------------------------- FILE CONVERSION AND SORTING ---------------------------"
        
        # Convert SeqAlto alignment to BAM and sort
        seq_bam_conversion_start=$(timer)
        
        echo "Converting SeqAlto alignment to BAM."
        samtools view -bS ${reads_prefix}_SeqAlto.sam > ${reads_prefix}_SeqAlto.bam
        
        echo "SeqAlto SAM-to-BAM conversion complete in $(timer ${seq_bam_conversion_start})." | tee -a ${log_file}
        
        # Convert BWA alignment to BAM and sort
        bwa_bam_conversion_start=$(timer)
        
        echo "Converting BWA alignment to BAM."
        samtools view -bS ${reads_prefix}_BWA.sam > ${reads_prefix}_BWA.bam
        
        echo "BWA SAM-to-BAM conversion complete in $(timer ${bwa_bam_conversion_start})." | tee -a ${log_file}

        echo -e "\n--------------------------- RUNNING ALIGNSTATS ---------------------------"
        
        alignstats_start=$(timer)
        
        # Run AlignDiff
        java -jar ${alignstats} \
            ${reads_prefix}_SeqAlto.bam ${reads_prefix}_BWA.bam \
            ${reads_prefix}_compare/ \
            --report ${reads_prefix}_compare.html \
            --temp-dir ${tmp}
        #    --pivot-report ${results_file} \
            
        echo "AlignStats analysis complete in $(timer ${alignstats_start})."
        echo "Alignment Difference report written to '${reads_prefix}_compare.html'."
    
        ;;

    * )
        # Do not perform runstats
        ;;

esac

echo -e "\n--------------------------- ALL ALIGNMENTS COMPLETE ---------------------------"
