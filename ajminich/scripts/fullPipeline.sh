#!/bin/bash -eu

# Full Pipeline Runner
# Written by AJ Minich (aj.minich@binatechnologies.com), February 2012
#
# Runs alignment and the full GATK pipeline on a given alignment file.
#

# Constants
PICARD=~/programs/picard/dist
GATK=~/programs/gatk/dist
seqalto=/home/ajminich/programs/seq
bwa=/home/ajminich/programs/bwa/bwa-0.6.1/bwa
alignstats=/home/ajminich/programs/alignstats
SNP=/mnt/scratch0/public/data/variants/dbSNP/dbsnp_132.hg19.vcf

tmp=/mnt/scratch0/ajminich
aligner_threads=16
pipeline_threads=1

if [[ $# -lt 4 ]]; then
    echo "Full Pipeline Runner v1.0"
    echo "Usage: fullPipeline <reference> <reads_prefix> <aligner> <pipeline>"
    echo ""
    echo "Reference:            location of the reference FASTA file, pre-indexed"
    echo "Reads prefix:         reads will be taken from <reads_prefix>_1.fq and <reads_prefix>_2.fq"
    echo "Aligner choices:      seqalto, bwa"
    echo "Pipeline choices:     gatk"
    echo ""
    echo "Notes:"
    echo "- the alignment is run with ${aligner_threads} thread(s)."
    echo "- the pipeline is run with ${pipeline_threads} thread(s)."
    exit
else
    reference=${1}
    reads=${2}
    aligner=${3}
    pipeline=${4}
    
    # Determine chromosome using the reference filename
    ref_file=$(basename $reference)
    chrom=${ref_file%.*}
    
    echo "Running ${aligner} alignment and ${pipeline} analysis on ${chrom} reads from '${reads}' with ${reference}."
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
        
        bwa_time=`echo $(timer ${bwa_start})`
        echo "BWA alignment complete in ${bwa_time}."
        
        ;;
        
    * )
        echo "Error: I don't know how to use '${aligner}'."
        exit
        ;;

esac

align_time=`echo $(timer ${align_start})`
echo "Alignment with ${aligner} complete in ${align_time}."

aligned_reads=${reads}\_${aligner}

echo -e "\n--------------------------- SAM-to-BAM CONVERSION AND SORTING ---------------------------"

conversion_start=$(timer)

echo "Converting SAM results file to BAM."
samtools view -bS ${aligned_reads}.sam > ${aligned_reads}.bam

conversion_time=`echo $(timer ${conversion_start})`
echo "SAM-to-BAM conversion complete in ${conversion_time}."

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

groups_time=`echo $(timer ${groups_start})`
echo "Add/Replace Groups complete in ${groups_time}."

aligned_reads=${aligned_reads}\_AG

sorting_start=$(timer)

java -Xms5g -Xmx5g -jar ${PICARD}/SortSam.jar \
    INPUT=${aligned_reads}.bam \
    VALIDATION_STRINGENCY=LENIENT \
    OUTPUT=${aligned_reads}\_sorted.bam \
    SORT_ORDER=coordinate

sorting_time=`echo $(timer ${sorting_start})`
echo "BAM sorting complete in ${sorting_time}."

aligned_reads=${aligned_reads}\_sorted

# Rebuild index
#java -Xms5g -Xmx5g -jar ${PICARD}/BuildBamIndex.jar INPUT=$d/$f.bam VALIDATION_STRINGENCY=LENIENT
#java -Xms15g -Xmx15g -jar $PICARD/ReorderSam.jar I=$d/$f.bam O=$d/$f\_sorted.bam  REFERENCE=$REF TMP_DIR=$TMP VALIDATION_STRINGENCY=LENIENT

echo -e ""
echo -e "Timing Results to this point:"
echo -e "${aligner} Alignment:      ${align_time}"
echo -e "SAM-to-BAM Conversion:     ${conversion_time}"
echo -e "Add/Remove Groups:         ${groups_time}"
echo -e "BAM Sorting:               ${sorting_time}"

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

mark_duplicates_time=`echo $(timer ${mark_duplicates_start})`
echo "Marking duplicates complete in ${mark_duplicates_time}."

echo -e "\n--------------------------- REALIGNMENT ---------------------------"

realigner_target_creator_start=$(timer)

echo "Determining intervals to re-align."
java -Xms5g -Xmx5g -jar ${GATK}/GenomeAnalysisTK.jar \
    -T RealignerTargetCreator \
    -I ${aligned_reads}\_marked.bam \
    -R ${reference} \
    -o ${aligned_reads}\_realign.intervals \
    -et NO_ET \
    -nt ${pipeline_threads} \
    -L ${chrom}

realigner_target_creator_time=`echo $(timer ${realigner_target_creator_start})`
echo "Realigner Target Creator complete in ${realigner_target_creator_time}."

indel_realigner_start=$(timer)

echo "Realigning targeted intervals."
java -Xms5g -Xmx5g -Djava.io.tmpdir=${tmp} -jar $GATK/GenomeAnalysisTK.jar \
    -T IndelRealigner \
    -I ${aligned_reads}\_marked.bam \
    -R ${reference} \
    -o ${aligned_reads}\_realign.bam \
    -targetIntervals ${aligned_reads}\_realign.intervals \
    -et NO_ET \
    -L ${chrom}

# Rebuild index
java -Xms5g -Xmx5g -jar $PICARD/BuildBamIndex.jar INPUT=${aligned_reads}\_realign.bam VALIDATION_STRINGENCY=LENIENT

indel_realigner_time=`echo $(timer ${indel_realigner_start})`
echo "Indel Realigner complete in ${indel_realigner_time}."

echo "Note: fix mate information executed by Indel Realigner."

echo -e ""
echo -e "Timing Results to this point:"
echo -e "${aligner} Alignment:      ${align_time}"
echo -e "SAM-to-BAM Conversion:     ${conversion_time}"
echo -e "Add/Remove Groups:         ${groups_time}"
echo -e "BAM Sorting:               ${sorting_time}"
echo -e "Mark Duplicates:           ${mark_duplicates_time}"
echo -e "Realigner Target Creator:  ${realigner_target_creator_time}"
echo -e "Indel Realigner:           ${indel_realigner_time}"

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
    -nt ${pipeline_threads} \
    -L ${chrom}

# Rebuild index
java -Xms5g -Xmx5g -jar $PICARD/BuildBamIndex.jar INPUT=${aligned_reads}\_realign.bam VALIDATION_STRINGENCY=LENIENT

count_covariates_time=`echo $(timer ${count_covariates_start})`
echo "Count Covariates complete in ${count_covariates_time}."

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
        -et NO_ET \
        -L ${chrom}
fi

# Rebuild index
java -Xms5g -Xmx5g -jar $PICARD/BuildBamIndex.jar INPUT=${aligned_reads}\_recalibrated.bam VALIDATION_STRINGENCY=LENIENT

recalibrate_table_time=`echo $(timer ${recalibrate_table_start})`
echo "Table Recalibration complete in ${recalibrate_table_time}."

echo -e "\n--------------------------- UNIFIED GENOTYPER SNP CALLING ---------------------------"

snp_calling_start=$(timer)

echo "Running GATK SNP caller."
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
    -nt ${pipeline_threads} \
    -L ${chrom}

snp_calling_time=`echo $(timer ${snp_calling_start})`
echo "SNP calling complete in ${snp_calling_time}."

indel_calling_start=$(timer)

echo "Running GATK indel caller."
java -Xms5g -Xmx5g -Djava.io.tmpdir=${tmp} -jar $GATK/GenomeAnalysisTK.jar \
    -T UnifiedGenotyper \
    -I ${aligned_reads}\_recalibrated.bam \
    -R ${reference} \
    -D ${SNP} \
    -o ${aligned_reads}\_indel.vcf \
    -et NO_ET \
    -A AlleleBalance \
    -A DepthOfCoverage \
    -A MappingQualityZero \
    -baq CALCULATE_AS_NECESSARY \
    -stand_call_conf 30.0 \
    -stand_emit_conf 10.0 \
    -glm INDEL \
    -nt ${pipeline_threads} \
    -L ${chrom}

indel_calling_time=`echo $(timer ${indel_calling_start})`
echo "Indel calling complete in ${indel_calling_time}."

echo -e "\n--------------------------- PIPELINE PROCESSING COMPLETE ---------------------------"

total_time=`echo $(timer ${start_time})`

echo -e "${aligner} Alignment:      ${align_time}"
echo -e "SAM-to-BAM Conversion:     ${conversion_time}"
echo -e "Add/Remove Groups:         ${groups_time}"
echo -e "BAM Sorting:               ${sorting_time}"
echo -e "Mark Duplicates:           ${mark_duplicates_time}"
echo -e "Realigner Target Creator:  ${realigner_target_creator_time}"
echo -e "Indel Realigner:           ${indel_realigner_time}"
echo -e "Counting Covariates:       ${count_covariates_time}"
echo -e "Table Recalibration:       ${recalibrate_table_time}"
echo -e "SNP Calling:               ${snp_calling_time}"
echo -e "Indel Calling:             ${indel_calling_time}"
echo -e "Total Running Time:        ${total_time}"