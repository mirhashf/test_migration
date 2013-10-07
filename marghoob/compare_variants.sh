#!/bin/bash

set -xe
export PERL5LIB=$HOME/vcftools_0.1.11/lib/perl5/site_perl
export JAVA_HOME=$HOME/lake/opt/jdk1.7.0_25/
export PATH=$HOME/vcftools_0.1.11/bin:$JAVA_HOME/bin:$PATH
export GATK_JAR=$HOME/lake/opt/gatk-2.6-5-gba531bd/GenomeAnalysisTK.jar

dir1=$1
dir2=$2
outdir=$3
reference=$4
outfile=$5

rm -f $outfile

DIR="$( cd "$( dirname "$0" )" && pwd )"

mkdir -pv $outdir
awk -v dir1=$dir1 '{print dir1"/vcfs/"$1".vcf.gz"}' $reference.fai > $outdir/files1
awk -v dir2=$dir2 '{print dir2"/vcfs/"$1".vcf.gz"}' $reference.fai > $outdir/files2

echo "Concatenating the chromosome VCFs"
vcf1=$outdir/all1.vcf.gz
vcf2=$outdir/all2.vcf.gz

(vcf-concat -f $outdir/files1 | bgzip > $vcf1; tabix -f $vcf1) &
(vcf-concat -f $outdir/files2 | bgzip > $vcf2; tabix -f $vcf2) &
wait

echo "Separating into indels and SNPs"
for vartype in SNP INDEL
do

(java -Xmx1g -Xms1g -jar $GATK_JAR -T SelectVariants -U LENIENT_VCF_PROCESSING -selectType $vartype -V $vcf1 -o $outdir/"$vartype"1.vcf -R $reference --excludeFiltered &>$outdir/"$vartype"1.log; bgzip -f $outdir/"$vartype"1.vcf; tabix -f $outdir/"$vartype"1.vcf.gz) &
(java -Xmx1g -Xms1g -jar $GATK_JAR -T SelectVariants -U LENIENT_VCF_PROCESSING -selectType $vartype -V $vcf2 -o $outdir/"$vartype"2.vcf -R $reference --excludeFiltered &>$outdir/"$vartype"2.log; bgzip -f $outdir/"$vartype"2.vcf; tabix -f $outdir/"$vartype"2.vcf.gz) &

done
wait

echo "Comparing the indel and SNP vcfs and getting stats"
for vartype in SNP INDEL
do
vcf-compare -a $outdir/"$vartype"1.vcf.gz $outdir/"$vartype"2.vcf.gz > $outdir/$vartype.vcf-compare.txt &
vcf-stats $outdir/$vartype"1.vcf.gz" -p $outdir/$vartype.1.stats &
vcf-stats $outdir/$vartype"2.vcf.gz" -p $outdir/$vartype.2.stats &
done

echo "Generating known and novel subsets of INDELS and SNPs"
for vartype in SNP INDEL
do

mkdir -p $outdir/$vartype

f1=$outdir/"$vartype"1.vcf.gz
f2=$outdir/"$vartype"2.vcf.gz

(vcf-isec $f1 $f2 | bgzip > $outdir/$vartype/common.vcf.gz; tabix -f $outdir/$vartype/common.vcf.gz) &
(vcf-isec -c $f1 $f2 | bgzip >  $outdir/$vartype/1.vcf.gz; tabix -f $outdir/$vartype/1.vcf.gz) &
(vcf-isec -c $f2 $f1 | bgzip > $outdir/$vartype/2.vcf.gz; tabix -f $outdir/$vartype/2.vcf.gz) &
wait

for subset in common 1 2
do

mkdir -p $outdir/$vartype/$subset
(cat <(gunzip -c $outdir/$vartype/$subset.vcf.gz|grep "^#") <(gunzip -c $outdir/$vartype/$subset.vcf.gz |grep -v "^#"|awk '{if ($3 != ".") print $0}') | bgzip > $outdir/$vartype/$subset/known.vcf.gz; tabix -f $outdir/$vartype/$subset/known.vcf.gz) &
(cat <(gunzip -c $outdir/$vartype/$subset.vcf.gz|grep "^#") <(gunzip -c $outdir/$vartype/$subset.vcf.gz |grep -v "^#"|awk '{if ($3 == ".") print $0}') | bgzip > $outdir/$vartype/$subset/novel.vcf.gz; tabix -f $outdir/$vartype/$subset/novel.vcf.gz) &
wait

for subsubset in known novel
do
vcf-stats $outdir/$vartype/$subset/$subsubset.vcf.gz -p $outdir/$vartype/$subset/$subsubset &
done

vcf-stats $outdir/$vartype/$subset.vcf.gz -p $outdir/$vartype/$subset.stats &

done

done

wait

echo "Getting the different counts"
$DIR/gen_tables.sh $outdir $dir1 $dir2 > $outfile
