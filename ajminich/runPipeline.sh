#!/bin/bash -eu

# Pipeline Runner
# Written by Jian Li
# Modified by AJ Minich (aj.minich@binatechnologies.com), January 2012
#
# Runs the full GATK pipeline on a given alignment file.
#
# To run the pipeline on the result of an alignment - for example,
# "alignment_result.sam" - under your home directory "/home/ajminich/", run:
#
#    runPipeline.sh /home/ajminich alignment_result
#
# Validation: http://www.broadinstitute.org/gsa/gatkdocs/release/org_broadinstitute_sting_gatk_walkers_varianteval_VariantEvalWalker.html

# d = home directory
d=$1

# f = SAM result file
f=$2

# Genome
chrom=chr1
#REF=/mnt/scratch0/public/data/annotation/hg19/ucsc.hg19.fasta     # run against full genome
REF=/mnt/scratch0/public/data/annotation/hg19/chr1.fa     # run against chr1
SNP=/mnt/scratch0/public/data/variants/dbSNP/dbsnp_132.hg19.vcf

# Constants
PICARD=~/programs/picard/dist
GATK=~/programs/gatk/dist
tmp=/mnt/scratch0/ajminich
num_threads=1

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

# Clean up current directory (this could be dangerous!)
#rm -f ${f}*.[bam,bai,metrics,intervals,vcf,idx]

start_time=$(timer)

echo -e "\n--------------------------- SAM-to-BAM CONVERSION ---------------------------"

file_conversion_start=$(timer)

#samtools view -bt $REF.fai $d/$f.sam > $d/$f.bam
#samtools view -bS $d/$f.sam > $d/$f.bam
#java -Xms15g -Xmx15g -jar $PICARD/MergeSamFiles.jar INPUT=$rg/1_RG.bam INPUT=$rg/2_RG.bam INPUT=$rg/3_RG.bam INPUT=$rg/4_RG.bam INPUT=$rg/5_RG.bam OUTPUT=$rg/bwa_merged.bam SORT_ORDER=coordinate VALIDATION_STRINGENCY=LENIENT  ASSUME_SORTED=true USE_THREADING=true

#java -Xms5g -Xmx5g -jar $PICARD/AddOrReplaceReadGroups.jar \
#    I=$d/$f.bam \
#    O=$d/$f\_AG.bam \
#    SORT_ORDER=coordinate \
#    RGID=SIMUG21 \
#    RGLB=READGROUPLIBRARY1 \
#    RGPL=illumina \
#    RGSM=Blah RGPU=hi \
#    VALIDATION_STRINGENCY=LENIENT  \
#    TMP_DIR=${tmp}

f=${f}\_AG

#java -Xms5g -Xmx5g -jar $PICARD/SortSam.jar \
#    INPUT=$d/$f.bam \
#    VALIDATION_STRINGENCY=LENIENT \
#    OUTPUT=$d/$f\_sorted.bam \
#    SORT_ORDER=coordinate

f=$f\_sorted

# Rebuild index
#java -Xms5g -Xmx5g -jar ${PICARD}/BuildBamIndex.jar INPUT=$d/$f.bam VALIDATION_STRINGENCY=LENIENT
#java -Xms15g -Xmx15g -jar $PICARD/ReorderSam.jar I=$d/$f.bam O=$d/$f\_sorted.bam  REFERENCE=$REF TMP_DIR=$TMP VALIDATION_STRINGENCY=LENIENT

file_conversion_time=`echo $(timer ${file_conversion_start})`
echo "SAM-to-BAM conversion and sorting complete in ${file_conversion_time}."

echo -e "\n--------------------------- DUPLICATE MARKING ---------------------------"

mark_duplicates_start=$(timer)

echo "Marking duplicates."
#java -Xms5g -Xmx5g -jar ${PICARD}/MarkDuplicates.jar \
#    TMP_DIR=${tmp} \
#    I=$d/$f.bam\
#    O=$d/$f\_marked.bam\
#    M=$d/$f.metrics \
#    VALIDATION_STRINGENCY=SILENT \
#    ASSUME_SORTED=true \
#    REMOVE_DUPLICATES=true 

# Rebuild index
java -Xms5g -Xmx5g -jar ${PICARD}/BuildBamIndex.jar INPUT=$d/$f\_marked.bam VALIDATION_STRINGENCY=LENIENT

mark_duplicates_time=`echo $(timer ${mark_duplicates_start})`
echo "Marking duplicates complete in ${mark_duplicates_time}."

echo -e "\n--------------------------- REALIGNMENT ---------------------------"

realigner_target_creator_start=$(timer)

echo "Determining intervals to re-align."
java -Xms5g -Xmx5g -jar $GATK/GenomeAnalysisTK.jar \
    -T RealignerTargetCreator \
    -I $d/$f\_marked.bam \
    -R $REF \
    -o $d/$f\_realign.intervals \
    -et NO_ET \
    -nt ${num_threads} \
    -L $chrom

realigner_target_creator_time=`echo $(timer ${realigner_target_creator_start})`
echo "Realigner Target Creator complete in ${realigner_target_creator_time}."

indel_realigner_start=$(timer)

echo "Realigning targeted intervals."
java -Xms5g -Xmx5g -Djava.io.tmpdir=${tmp} -jar $GATK/GenomeAnalysisTK.jar \
    -T IndelRealigner \
    -I $d/$f\_marked.bam \
    -R $REF \
    -o $d/$f\_realign.bam \
    -targetIntervals $d/$f\_realign.intervals \
    -et NO_ET \
    -L $chrom
#        -nt ${num_threads} \       # Does not currently run in parallel mode

indel_realigner_time=`echo $(timer ${indel_realigner_start})`
echo "Indel Realigner complete in ${indel_realigner_time}."

# Rebuild index
java -Xms5g -Xmx5g -jar $PICARD/BuildBamIndex.jar INPUT=$d/$f\_realign.bam VALIDATION_STRINGENCY=LENIENT

mate_information_start=$(timer)

echo "Fixing mate pairs and order of realigned reads."
java -Xms5g -Xmx5g -jar $PICARD/FixMateInformation.jar \
    TMP_DIR=${tmp} \
    INPUT=$d/$f\_realign.bam \
    VALIDATION_STRINGENCY=SILENT \
    SORT_ORDER=coordinate

# Rebuild index
java -Xms5g -Xmx5g -jar $PICARD/BuildBamIndex.jar INPUT=$d/$f\_realign.bam VALIDATION_STRINGENCY=LENIENT

mate_information_time=`echo $(timer ${mate_information_start})`
echo "Fix Mate Information complete in ${mate_information_time}."

echo -e "\n--------------------------- BASE QUALITY RECALIBRATION ---------------------------"

count_covariates_start=$(timer)

echo "Counting covariates."
java -Xms5g -Xmx5g -jar $GATK/GenomeAnalysisTK.jar \
    -T CountCovariates \
    -cov ReadGroupCovariate \
    -cov QualityScoreCovariate \
    -cov CycleCovariate \
    -cov DinucCovariate \
    -R $REF \
    -knownSites ${SNP}\
    -I $d/$f\_realign.bam \
    -recalFile $d/$f\_realign.bam.csv \
    -et NO_ET \
    -nt ${num_threads} \
    -L $chrom

count_covariates_time=`echo $(timer ${count_covariates_start})`
echo "Count Covariates complete in ${count_covariates_time}."

# Rebuild index
java -Xms5g -Xmx5g -jar $PICARD/BuildBamIndex.jar INPUT=$d/$f\_realign.bam VALIDATION_STRINGENCY=LENIENT

recalibrate_table_start=$(timer)

echo "Recalibrating base quality table."
if [ "`grep -v '#' $d/$f\_realign.bam.csv | grep -v "EOF" | wc -l`" = "1" ]
then
    cp $d/$f\_realign.bam $d/$f\_recalibrated.bam
else
    java -Xms5g -Xmx5g -jar ${GATK}/GenomeAnalysisTK.jar \
        -R $REF \
        -I $d/$f\_realign.bam  \
        -o $d/$f\_recalibrated.bam \
        -T TableRecalibration \
        -baq RECALCULATE \
        --doNotWriteOriginalQuals \
        -recalFile $d/$f\_realign.bam.csv \
        -et NO_ET \
        -L $chrom
        #-nt ${num_threads} \
fi

# Rebuild index
java -Xms5g -Xmx5g -jar $PICARD/BuildBamIndex.jar INPUT=$d/$f\_recalibrated.bam VALIDATION_STRINGENCY=LENIENT

recalibrate_table_time=`echo $(timer ${recalibrate_table_start})`
echo "Table Recalibration complete in ${recalibrate_table_time}."

echo -e "\n--------------------------- UNIFIED GENOTYPER SNP CALLING ---------------------------"

snp_calling_start=$(timer)

echo "Running GATK SNP caller."
java -Xms5g -Xmx5g -Djava.io.tmpdir=${tmp} -jar $GATK/GenomeAnalysisTK.jar \
    -T UnifiedGenotyper \
    -I $d/$f\_recalibrated.bam \
    -R $REF \
    -D $SNP \
    -o $d/$f.vcf  \
    -et NO_ET \
    -dcov 1000 \
    -A AlleleBalance \
    -A DepthOfCoverage \
    -A MappingQualityZero \
    -baq CALCULATE_AS_NECESSARY \
    -stand_call_conf 30.0 \
    -stand_emit_conf 10.0 \
    -nt ${num_threads} \
    -L $chrom

snp_calling_time=`echo $(timer ${snp_calling_start})`
echo "SNP calling complete in ${snp_calling_time}."

indel_calling_start=$(timer)

echo "Running GATK indel caller."
java -Xms5g -Xmx5g -Djava.io.tmpdir=${tmp} -jar $GATK/GenomeAnalysisTK.jar \
    -T UnifiedGenotyper \
    -I $d/$f\_recalibrated.bam \
    -R $REF \
    -D $SNP \
    -o $d/$f\_indel.vcf \
    -et NO_ET \
    -A AlleleBalance \
    -A DepthOfCoverage \
    -A MappingQualityZero \
    -baq CALCULATE_AS_NECESSARY \
    -stand_call_conf 30.0 \
    -stand_emit_conf 10.0 \
    -glm INDEL \
    -nt ${num_threads} \
    -L $chrom

indel_calling_time=`echo $(timer ${indel_calling_start})`
echo "Indel calling complete in ${indel_calling_time}."

echo -e "\n--------------------------- PIPELINE PROCESSING COMPLETE ---------------------------"

total_time=`echo $(timer ${start_time})`

echo -e "File Conversion:           ${file_conversion_time}"
echo -e "Mark Duplicates:           ${mark_duplicates_time}"
echo -e "Realigner Target Creator:  ${realigner_target_creator_time}"
echo -e "Indel Realigner:           ${indel_realigner_time}"
echo -e "Fix Mate Information:      ${mate_information_time}"
echo -e "Counting Covariates:       ${count_covariates_time}"
echo -e "Table Recalibration:       ${recalibrate_table_time}"
echo -e "SNP Calling:               ${snp_calling_time}"
echo -e "Indel Calling:             ${indel_calling_time}"
echo -e "Total Running Time:        ${total_time}"