#!/bin/bash

set -xe 
export PERL5LIB=$HOME/vcftools_0.1.11/lib/perl5/site_perl
export JAVA_HOME=$HOME/lake/opt/jdk1.7.0_25/
export PATH=$HOME/vcftools_0.1.11/bin:$JAVA_HOME/bin:$PATH
export GATK_JAR=$HOME/lake/opt/gatk-2.6-5-gba531bd/GenomeAnalysisTK.jar

jobdir=$1
tmpdir=$2
reference=$3
outfile=$4

rm -f $outfile

DIR="$( cd "$( dirname "$0" )" && pwd )"

mkdir -pv $tmpdir
awk -v jobdir=$jobdir '{print jobdir"/vcfs/"$1".vcf.gz"}' $reference.fai > $tmpdir/files

echo "Concatenating the chromosome VCFs"
vcf=$tmpdir/all.vcf.gz

(vcf-concat -f $tmpdir/files | bgzip > $vcf; tabix -f $vcf)

declare -A exclude_filter
exclude_filter["PASS"]="--excludeFiltered"
exclude_filter["ALL"]=""
exclude_filter["NONPASS"]=""

for filter in ALL PASS NONPASS
do

echo "Separating into indels and SNPs"
if [ "$filter" == "NONPASS" ]
then

cat <(gunzip -c $vcf|grep "^#") <(gunzip -c $vcf|grep -v "^#"|awk '{if ($7 != "PASS") print $0}') |bgzip > $tmpdir/nopass.vcf.gz
tabix -f $tmpdir/nopass.vcf.gz

for vartype in SNP INDEL
do
(java -jar -Xmx1g -Xms1g $GATK_JAR -T SelectVariants -U LENIENT_VCF_PROCESSING -selectType $vartype -V $tmpdir/nopass.vcf.gz -o $tmpdir/"$vartype".vcf -R $reference ${exclude_filter[$filter]} &>$tmpdir/"$vartype".log; bgzip -f $tmpdir/"$vartype".vcf; tabix -f $tmpdir/"$vartype".vcf.gz) &
done

else
for vartype in SNP INDEL
do
(java -jar -Xmx1g -Xms1g $GATK_JAR -T SelectVariants -U LENIENT_VCF_PROCESSING -selectType $vartype -V $vcf -o $tmpdir/"$vartype".vcf -R $reference ${exclude_filter[$filter]} &>$tmpdir/"$vartype".log; bgzip -f $tmpdir/"$vartype".vcf; tabix -f $tmpdir/"$vartype".vcf.gz) &
done
fi
wait

echo "Getting stats on SNPs and indels"
for vartype in SNP INDEL
do
vcf-stats $tmpdir/$vartype".vcf.gz" -p $tmpdir/$vartype &
done

echo "Generating known and novel subsets of INDELS and SNPs"
for vartype in SNP INDEL
do

mkdir -p $tmpdir/$vartype

(cat <(gunzip -c $tmpdir/$vartype.vcf.gz|grep "^#") <(gunzip -c $tmpdir/$vartype.vcf.gz |grep -v "^#"|awk '{if ($3 != ".") print $0}') | bgzip > $tmpdir/$vartype/known.vcf.gz; tabix -f $tmpdir/$vartype/known.vcf.gz) &
(cat <(gunzip -c $tmpdir/$vartype.vcf.gz|grep "^#") <(gunzip -c $tmpdir/$vartype.vcf.gz |grep -v "^#"|awk '{if ($3 == ".") print $0}') | bgzip > $tmpdir/$vartype/novel.vcf.gz; tabix -f $tmpdir/$vartype/novel.vcf.gz) &

done
wait

for vartype in SNP INDEL
do
for subset in known novel
do
vcf-stats $tmpdir/$vartype/$subset.vcf.gz -p $tmpdir/$vartype/$subset &
done
done

wait

echo "Getting the different counts"
$DIR/gen_single_tables.sh $tmpdir "$jobdir-$filter" >> $outfile
echo "" >> $outfile

done

#rm -rf $tmpdir
