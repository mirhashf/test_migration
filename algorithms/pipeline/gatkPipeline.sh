#!/bin/bash -eu


echo "*** GATK Pipeline ***"

SNP=/mnt/scratch0/public/data/variants/dbSNP/dbsnp_132.hg19.vcf
PICARD=~jianl/downloads/picard/dist
GATK=~jianl/downloads/gatk/dist
REF=/mnt/scratch0/public/data/annotation/hg19/chr21.fa

# To run the pipeline on the result of alignment, for example, the "result.sam"
# you got from running seqAlto under your home directory "/home/ajminich/", run:
# gatkPipline.sh /home/ajminich result
# f = SAM result file
# d = home directory
f=$2
d=$1
chrom=chr21

mkdir  -p $d/tmp
TMP=$d/tmp

#samtools view -bt $REF.fai $d/$f.sam > $d/$f.bam
samtools view -bS $d/$f.sam > $d/$f.bam
#java -Xms15g -Xmx15g -jar $PICARD/MergeSamFiles.jar INPUT=$rg/1_RG.bam INPUT=$rg/2_RG.bam INPUT=$rg/3_RG.bam INPUT=$rg/4_RG.bam INPUT=$rg/5_RG.bam OUTPUT=$rg/bwa_merged.bam SORT_ORDER=coordinate VALIDATION_STRINGENCY=LENIENT  ASSUME_SORTED=true USE_THREADING=true

java -Xms5g -Xmx5g -jar $PICARD/AddOrReplaceReadGroups.jar I=$d/$f.bam O=$d/$f\_AG.bam SORT_ORDER=coordinate RGID=SIMUG21 RGLB=READGROUPLIBRARY1 RGPL=illumina RGSM=Blah RGPU=hi VALIDATION_STRINGENCY=LENIENT  TMP_DIR=$TMP
f=$f\_AG

java -Xms5g -Xmx5g -jar $PICARD/SortSam.jar INPUT=$d/$f.bam VALIDATION_STRINGENCY=LENIENT OUTPUT=$d/$f\_sorted.bam SORT_ORDER=coordinate

f=$f\_sorted

java -Xms5g -Xmx5g -jar $PICARD/BuildBamIndex.jar INPUT=$d/$f.bam VALIDATION_STRINGENCY=LENIENT
#java -Xms15g -Xmx15g -jar $PICARD/ReorderSam.jar I=$d/$f.bam O=$d/$f\_sorted.bam  REFERENCE=$REF TMP_DIR=$TMP VALIDATION_STRINGENCY=LENIENT

echo ">>> Marking duplicates"

java -Xms5g -Xmx5g -jar $PICARD/MarkDuplicates.jar \
        TMP_DIR=$TMP \
        I=$d/$f.bam\
        O=$d/$f\_marked.bam\
        M=$d/$f.metrics \
        VALIDATION_STRINGENCY=SILENT \
        ASSUME_SORTED=true \
        REMOVE_DUPLICATES=true 


echo "*** Finished removing duplicates ***"

java -Xms5g -Xmx5g -jar $PICARD/BuildBamIndex.jar INPUT=$d/$f\_marked.bam VALIDATION_STRINGENCY=LENIENT

echo ">>> Determining (small) suspicious intervals which are likely in need of realignment"
java -Xms5g -Xmx5g -jar $GATK/GenomeAnalysisTK.jar \
        -T RealignerTargetCreator \
        -I $d/$f\_marked.bam \
        -R $REF \
        -o $d/$f\_realign.intervals \
        -et NO_ET \
 -L $chrom

echo ">>> Running the realigner over the targeted intervals"
java -Xms5g -Xmx5g -Djava.io.tmpdir=$TMP -jar $GATK/GenomeAnalysisTK.jar \
        -T IndelRealigner \
        -I $d/$f\_marked.bam \
        -R $REF \
        -o $d/$f\_realign.bam \
        -targetIntervals $d/$f\_realign.intervals \
        -et NO_ET \
 -L $chrom

java -Xms5g -Xmx5g -jar $PICARD/BuildBamIndex.jar INPUT=$d/$f\_realign.bam VALIDATION_STRINGENCY=LENIENT


echo ">>> Fixing the mate pairs and order of the realigned reads"
java -Xms5g -Xmx5g -jar $PICARD/FixMateInformation.jar \
        TMP_DIR=$TMP \
        INPUT=$d/$f\_realign.bam \
        VALIDATION_STRINGENCY=SILENT \
        SORT_ORDER=coordinate
echo "*** Finished realigning targeted regions ***" 

java -Xms5g -Xmx5g -jar $PICARD/BuildBamIndex.jar INPUT=$d/$f\_realign.bam VALIDATION_STRINGENCY=LENIENT

echo ">>> Counting covariates"
java -Xms5g -Xmx5g -jar $GATK/GenomeAnalysisTK.jar \
        -T CountCovariates \
        -cov ReadGroupCovariate \
        -cov QualityScoreCovariate \
        -cov CycleCovariate \
        -cov DinucCovariate \
        -R $REF \
    -knownSites $SNP\
        -I $d/$f\_realign.bam \
        -recalFile $d/$f\_realign.bam.csv \
        -et NO_ET \
        -nt 8 \
    -L $chrom


java -Xms5g -Xmx5g -jar $PICARD/BuildBamIndex.jar INPUT=$d/$f\_realign.bam VALIDATION_STRINGENCY=LENIENT

echo ">> Table recalibration"
if [ "`grep -v '#' $d/$f\_realign.bam.csv | grep -v "EOF" | wc -l`" = "1" ]
then
        cp $d/$f\_realign.bam $d/$f\_recalibrated.bam
else
        java -Xms5g -Xmx5g -jar $GATK/GenomeAnalysisTK.jar \
        -R $REF \
        -I $d/$f\_realign.bam  \
        -o $d/$f\_recalibrated.bam \
        -T TableRecalibration \
        -baq RECALCULATE \
        --doNotWriteOriginalQuals \
        -recalFile $d/$f\_realign.bam.csv \
        -et NO_ET \
	    -L $chrom

fi

echo "*** Finished recalibrating base quality ***"

java -Xms5g -Xmx5g -jar $PICARD/BuildBamIndex.jar INPUT=$d/$f\_recalibrated.bam VALIDATION_STRINGENCY=LENIENT

echo ">>> Running the unified genotyper for SNP calling"
java -Xms5g -Xmx5g -Djava.io.tmpdir=$TMP -jar $GATK/GenomeAnalysisTK.jar \
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
    -nt 2 \
    -L $chrom

echo "*** Finished SNP Analysis using the GATK Unified Genotyper ***"

java -Xms5g -Xmx5g -Djava.io.tmpdir=$TMP -jar $GATK/GenomeAnalysisTK.jar \
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
 -L $chrom 

