#!/bin/bash -eu

# Full Pipeline Runner
# Written by AJ Minich (aj.minich@binatechnologies.com), February 2012
#
# Runs alignment and a full GATK pipeline on a given alignment file.
# Logs output to <chrom>_<aligner>_<pipeline>.log
#

VERSION=1.02

# Constants
PICARD=/home/ajminich/programs/picard/dist
GATK=/home/ajminich/programs/gatk/dist
GATK_FAST=/home/ajminich/programs/third_party/software/gatk/dist/
seqalto=/home/ajminich/programs/seq
bwa=/home/ajminich/programs/bwa/bwa-0.6.1/bwa
SNP=/mnt/scratch0/public/data/variants/dbSNP/dbsnp_132.hg19.vcf

tmp=/mnt/scratch0/ajminich
aligner_threads=1
pipeline_threads=1

if [[ $# -lt 4 ]]; then
    echo "Full Pipeline Runner v${VERSION}"
    echo "Usage: pipeline.sh <reference> <reads_prefix> <aligner> <pipeline>"
    echo ""
    echo "Reference:            location of the reference FASTA file, pre-indexed"
    echo "Reads prefix:         reads will be taken from <reads_prefix>_1.fq and <reads_prefix>_2.fq"
    echo "Aligner choices:      seqalto, bwa"
    echo "Pipeline choices:     gatk, gatk-fast"
    echo ""
    echo "Notes:"
    echo "- Alignment is run with ${aligner_threads} thread(s)."
    echo "- Pipeline is run with ${pipeline_threads} thread(s)."
    echo "- Fix mate information is not executed separately, as it is performed by Indel Realigner."
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
    echo "Aligner: ${aligner} with ${aligner_threads} thread(s)" | tee -a ${log_file}
    echo "Pipeline: ${pipeline} with ${pipeline_threads} thread(s)" | tee -a ${log_file}
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

# Run the alignment
align_start=$(timer)
case "${aligner}" in

    "seqalto" | "seq" | "SeqAlto" | "Seqalto" )
    
        aligner="SeqAlto"
    
        echo -e "\n--------------------------- ALIGNMENT WITH SEQALTO ---------------------------"

        #index_load_start=$(timer)
        #${seqalto} load_index ${reference}\_22.midx
        #echo "SeqAlto Index Loading: $(timer ${index_load_start})" | tee -a ${log_file}

        ${seqalto} align \
            --idx ${reference}\_22.midx \
            -p ${aligner_threads} \
            -1 ${reads}\_1.fq \
            -2 ${reads}\_2.fq \
            > ${reads}\_${aligner}.sam
            
        ;;
        
    "bwa" | "BWA" )
    
        echo -e "\n--------------------------- ALIGNMENT WITH BWA ---------------------------"

        aligner="BWA"

        time ${bwa} aln \
            ${reference} \
            ${reads}\_1.fq \
            -t ${aligner_threads} \
            > ${reads}\_1.sai
        time ${bwa} aln \
            ${reference} \
            ${reads}\_2.fq \
            -t ${aligner_threads} \
            > ${reads}\_2.sai
        
        echo "Performing BWA paired-end merging."
        time ${bwa} sampe \
            -P ${reference} \
            ${reads}\_1.sai \
            ${reads}\_2.sai \
            ${reads}\_1.fq \
            ${reads}\_2.fq \
            > ${reads}\_${aligner}.sam
        
        ;;
        
    * )
        echo "Error: I don't know how to use '${aligner}'."
        exit
        ;;

esac

echo "${aligner} Alignment: $(timer ${align_start})" | tee -a ${log_file}

aligned_reads=${reads}\_${aligner}

echo -e "\n--------------------------- SAM-to-BAM CONVERSION AND SORTING ---------------------------"

conversion_start=$(timer)

samtools view -bS ${aligned_reads}.sam > ${aligned_reads}.bam

echo "SAM-to-BAM Conversion: $(timer ${conversion_start})" | tee -a ${log_file}

#java -Xms15g -Xmx15g -jar $PICARD/MergeSamFiles.jar INPUT=$rg/1_RG.bam INPUT=$rg/2_RG.bam INPUT=$rg/3_RG.bam INPUT=$rg/4_RG.bam INPUT=$rg/5_RG.bam OUTPUT=$rg/bwa_merged.bam SORT_ORDER=coordinate VALIDATION_STRINGENCY=LENIENT  ASSUME_SORTED=true USE_THREADING=true

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

# Rebuild index
#java -Xms5g -Xmx5g -jar ${PICARD}/BuildBamIndex.jar INPUT=$d/$f.bam VALIDATION_STRINGENCY=LENIENT
#java -Xms15g -Xmx15g -jar $PICARD/ReorderSam.jar I=$d/$f.bam O=$d/$f\_sorted.bam  REFERENCE=$REF TMP_DIR=$TMP VALIDATION_STRINGENCY=LENIENT

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

echo "Total Running Time: $(timer ${start_time})"
echo "Log written to '${log_file}'."

