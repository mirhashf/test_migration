#!/bin/bash -eu

# Pipeline Runner
# Written by AJ Minich (aj.minich@binatechnologies.com), February 2012
#
# Runs GATK pipeline on a given alignment file.
# Logs output to <alignment_file>_<pipeline>.log
#

VERSION=1.03

# Constants
PICARD=/home/ajminich/programs/picard/dist
GATK=/home/ajminich/programs/gatk/dist
GATK_FAST=/home/ajminich/programs/gatk-fast/software/gatk/dist/
SNP=/mnt/scratch0/public/data/variants/dbSNP/dbsnp_132.hg19.vcf

tmp=/tmp
pipeline_threads=16

if [[ $# -lt 4 ]]; then
    echo "Pipeline Runner v${VERSION}"
    echo "Usage: pipeline.sh <reference> <alignment_file>.sam <pipeline> <chroms>"
    echo ""
    echo "Reference:            location of the reference FASTA file, pre-indexed"
    echo "Alignment File:       alignment file to process (in .sam format)"
    echo "Pipeline choices:     gatk, gatk-fast"
    echo "Chromsosomes:         comma-separated list of chromosomes (chr15,chr16,...)"
    echo ""
    echo "Notes:"
    echo "- Pipeline is run with ${pipeline_threads} thread(s)."
    echo "- Fix mate information is not executed separately, as it is performed by Indel Realigner."
    exit
else
    reference=${1}
    alignment_file=${2}
    pipeline=${3}
    chrom=${4}
    
    # Initialize the log file
    log_file="${alignment_file}_${pipeline}.log"
    
    echo "Pipeline Runner v${VERSION}" | tee ${log_file}
    echo "Pipeline: ${pipeline} with ${pipeline_threads} thread(s)" | tee -a ${log_file}
    echo "Reference: ${reference}" | tee -a ${log_file}
    echo "Chromosomes: ${chrom}" | tee -a ${log_file}
    echo "Alignment File: ${alignment_file}" | tee -a ${log_file}
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

aligned_reads=${alignment_file}

echo -e "\n--------------------------- SAM-to-BAM CONVERSION AND SORTING ---------------------------"

#conversion_start=$(timer)

#samtools view -bS ${aligned_reads}.sam > ${aligned_reads}.bam

#echo "SAM-to-BAM Conversion: $(timer ${conversion_start})" | tee -a ${log_file}

groups_start=$(timer)

java -Xms5g -Xmx5g -jar ${PICARD}/AddOrReplaceReadGroups.jar \
    I=${aligned_reads}.bam \
    O=${aligned_reads}_AG.bam \
    SORT_ORDER=coordinate \
    RGID=SIMUG21 \
    RGLB=READGROUPLIBRARY1 \
    RGPL=illumina \
    RGSM=Blah RGPU=hi \
    VALIDATION_STRINGENCY=LENIENT  \
    TMP_DIR=${tmp}

echo "Add/Replace Groups: $(timer ${groups_start})" | tee -a ${log_file}

aligned_reads=${aligned_reads}\_AG

sorting_start=$(timer)

java -Xms5g -Xmx5g -jar ${PICARD}/SortSam.jar \
    INPUT=${aligned_reads}.bam \
    VALIDATION_STRINGENCY=LENIENT \
    OUTPUT=${aligned_reads}\_sorted.bam \
    SORT_ORDER=coordinate

echo "BAM Sorting: $(timer ${sorting_start})" | tee -a ${log_file}

aligned_reads=${aligned_reads}\_sorted

echo -e "\n--------------------------- DUPLICATE MARKING ---------------------------"

mark_duplicates_start=$(timer)

java -Xms5g -Xmx5g -jar ${PICARD}/MarkDuplicates.jar \
    TMP_DIR=${tmp} \
    I=${aligned_reads}.bam\
    O=${aligned_reads}\_marked.bam\
    M=${aligned_reads}.metrics \
    VALIDATION_STRINGENCY=SILENT \
    ASSUME_SORTED=true \
    REMOVE_DUPLICATES=true 

# Rebuild index
java -Xms5g -Xmx5g -jar ${PICARD}/BuildBamIndex.jar INPUT=${aligned_reads}\_marked.bam VALIDATION_STRINGENCY=LENIENT

echo "Mark Duplicates: $(timer ${mark_duplicates_start})" | tee -a ${log_file}

case "${pipeline}" in 

    "gatk-fast" | "fast-gatk" | "realigner-fast" | "fast" | "fastrealign" )
    
        pipeline="Fast GATK"
        
        echo -e "\n--------------------------- FAST GATK ---------------------------"
        
        fast_gatk_start=$(timer)
        
        java -Xms5g -Xmx5g -jar ${GATK_FAST}/GenomeAnalysisTK.jar \
            -T FastRealign \
            -I ${aligned_reads}\_marked.bam \
            -R ${reference} \
            -standard \
            -knownSites ${SNP} \
            -o ${aligned_reads}\_realign.bam \
            -recalFile ${aligned_reads}\_realign.bam.csv
            #--suppress_recalibrate
            #--suppress_write_bam
        
        echo "Fast GATK: $(timer ${fast_gatk_start})" | tee -a ${log_file}
        
        ;;

    * )     # everything else

        pipeline="Normal GATK"
        
        echo -e "\n--------------------------- REALIGNMENT ---------------------------"
        
        realigner_target_creator_start=$(timer)
        
        echo "Determining intervals to re-align."
        java -Xms5g -Xmx5g -jar ${GATK}/GenomeAnalysisTK.jar \
            -T RealignerTargetCreator \
            -I ${aligned_reads}\_marked.bam \
            -R ${reference} \
            -o ${aligned_reads}\_realign.intervals \
            -et NO_ET \
            -nt ${pipeline_threads}
        
        echo "Realigner Target Creator: $(timer ${realigner_target_creator_start})" | tee -a ${log_file}
        
        indel_realigner_start=$(timer)
        
        echo "Realigning targeted intervals."
        java -Xms5g -Xmx5g -Djava.io.tmpdir=${tmp} -jar $GATK/GenomeAnalysisTK.jar \
            -T IndelRealigner \
            -I ${aligned_reads}\_marked.bam \
            -R ${reference} \
            -o ${aligned_reads}\_realign.bam \
            -targetIntervals ${aligned_reads}\_realign.intervals \
            -et NO_ET
        
        # Rebuild index
        java -Xms5g -Xmx5g -jar $PICARD/BuildBamIndex.jar INPUT=${aligned_reads}\_realign.bam VALIDATION_STRINGENCY=LENIENT
        
        echo "Indel Realigner: $(timer ${indel_realigner_start})" | tee -a ${log_file}
        
        echo -e "\n--------------------------- BASE QUALITY RECALIBRATION ---------------------------"
        
        count_covariates_start=$(timer)
        
        echo "Counting covariates."
        java -Xms5g -Xmx5g -jar $GATK/GenomeAnalysisTK.jar \
            -T CountCovariates \
            -cov ReadGroupCovariate \
            -cov QualityScoreCovariate \
            -cov CycleCovariate \
            -cov DinucCovariate \
            -R ${reference} \
            -knownSites ${SNP} \
            -I ${aligned_reads}\_realign.bam \
            -recalFile ${aligned_reads}\_realign.bam.csv \
            -et NO_ET \
            -nt ${pipeline_threads}
        
        # Rebuild index
        java -Xms5g -Xmx5g -jar $PICARD/BuildBamIndex.jar INPUT=${aligned_reads}\_realign.bam VALIDATION_STRINGENCY=LENIENT
        
        echo "Counting Covariates: $(timer ${count_covariates_start})" | tee -a ${log_file}
        
        ;;

esac    # pipeline choice

echo -e "\n--------------------------- TABLE RECALIBRATION ---------------------------"

recalibrate_table_start=$(timer)

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

echo "Total Running Time: $(timer ${start_time})" | tee -a ${log_file}
echo "Log written to '${log_file}'."

