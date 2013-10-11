#!/bin/bash

set -xe 
jobdir=$1
tmpdir=$2
outfile=$3
gatk_ver=$4

declare -A gatk_paths
gatk_paths["1.6"]="$HOME/lake/opt/gatk-1.6"
gatk_paths["2.3-9"]="$HOME/lake/opt/gatk-2.3-9-gd785397"
gatk_paths["2.5-2"]="$HOME/lake/opt/CancerAnalysisPackage-2013.2-18-g8207e53"
gatk_paths["2.6-5"]="$HOME/lake/opt/gatk-2.6-5-gba531bd"
gatk_paths["2.7-2"]="$HOME/lake/opt/gatk-2.7-2-g6bda569"

GATK_PATH=${gatk_paths[$gatk_ver]}

export PERL5LIB=$HOME/vcftools_0.1.11/perl
export JAVA_HOME=$HOME/lake/opt/jdk1.7.0_25/
export PATH=$HOME/vcftools_0.1.11/bin:$JAVA_HOME/bin:$PATH
export GATK_JAR=$GATK_PATH/GenomeAnalysisTK.jar
export ANNOVAR=$HOME/lake/users/marghoob/annovar
export SNPSIFT=$HOME/lake/users/marghoob/snpEff/SnpSift.jar
export NISTVCF=$HOME/lake/users/marghoob/NIST/NISThighConf.hg19.vcf.gz
dbsnp=$HOME/lake/users/marghoob/GATK-bundle-hg19/dbsnp_137.hg19.vcf
reference=$HOME/lake/users/marghoob/GATK-bundle-hg19/ucsc.hg19.fa

ANNOTATE=true

rm -f $outfile

DIR="$( cd "$( dirname "$0" )" && pwd )"

mkdir -pv $tmpdir
awk -v jobdir=$jobdir '{print jobdir"/vcfs/"$1".vcf.gz"}' $reference.fai > $tmpdir/files

echo "Concatenating the chromosome VCFs"
vcf=$tmpdir/all.vcf.gz

if [ "$ANNOTATE" == "true" ]
then
echo "Concatenating and clearing ID field"
(vcf-concat -f $tmpdir/files |awk 'BEGIN{OFS="\t"} /^#/{print $0} !/^#/{printf("%s\t%s\t.\t",$1,$2); for(i=4;i<=NF;++i)printf("%s\t", $i);printf("\n")}'|bgzip > $tmpdir/all.pre_annotated.vcf.gz; tabix -f $tmpdir/all.pre_annotated.vcf.gz)
else
(vcf-concat -f $tmpdir/files | bgzip > $tmpdir/all.pre_annotated.vcf.gz; tabix -f $tmpdir/all.pre_annotated.vcf.gz)
fi

if [ "$ANNOTATE" == "true" ]
then
echo "Annotating variants"
(java -Xmx1g -Xms1g -jar $GATK_JAR -T VariantAnnotator -nt 8 -U LENIENT_VCF_PROCESSING -R $reference -D $dbsnp --variant $tmpdir/all.pre_annotated.vcf.gz --out $tmpdir/all.vcf &>$tmpdir/annotation.log; bgzip -f $tmpdir/all.vcf; tabix -f $tmpdir/all.vcf.gz)
else
mv $tmpdir/all.pre_annotated.vcf.gz $tmpdir/all.vcf.gz; mv $tmpdir/all.pre_annotated.vcf.gz.tbi $tmpdir/all.vcf.gz.tbi;
fi

declare -A exclude_filter
exclude_filter["PASS"]="--excludeFiltered"
exclude_filter["ALL"]=""
exclude_filter["NONPASS"]=""

for filter in PASS #NONPASS
do

echo "Separating into indels and SNPs"
if [ "$filter" == "NONPASS" ]
then

cat <(gunzip -c $vcf|grep "^#") <(gunzip -c $vcf|grep -v "^#"|awk '{if ($7 != "PASS") print $0}') |bgzip > $tmpdir/nopass.vcf.gz
tabix -f $tmpdir/nopass.vcf.gz

for vartype in SNP INDEL
do
(java -Xmx1g -Xms1g -jar $GATK_JAR -T SelectVariants -U LENIENT_VCF_PROCESSING -selectType $vartype -V $tmpdir/nopass.vcf.gz -o $tmpdir/"$vartype".vcf -R $reference ${exclude_filter[$filter]} &>$tmpdir/"$vartype".log; bgzip -f $tmpdir/"$vartype".vcf; tabix -f $tmpdir/"$vartype".vcf.gz) &
done

else
for vartype in SNP INDEL
do
(java -Xmx1g -Xms1g -jar $GATK_JAR -T SelectVariants -U LENIENT_VCF_PROCESSING -selectType $vartype -V $vcf -o $tmpdir/"$vartype".vcf -R $reference ${exclude_filter[$filter]} &>$tmpdir/"$vartype".log; bgzip -f $tmpdir/"$vartype".vcf; tabix -f $tmpdir/"$vartype".vcf.gz) &
done
fi
wait

if [ -n "$ANNOVAR" ]
then
echo "Converting to annovar format"
for vartype in SNP INDEL
do
$ANNOVAR/convert2annovar.pl -format vcf4 <(gunzip -c $tmpdir/$vartype.vcf.gz) -includeinfo -outfile $tmpdir/$vartype.avinput 2>$tmpdir/convert2annovar.log &
done
wait

echo "Annotating variants using annovar"
for vartype in SNP INDEL
do
$ANNOVAR/annotate_variation.pl --filter --dbtype snp137 -buildver hg19 $tmpdir/$vartype.avinput $ANNOVAR/humandb/ 2>$tmpdir/annotate_variation.log &
done
wait

for vartype in SNP INDEL
do
awk -v dbsnp=$dbsnp 'BEGIN{OFS=""; print "#!/bin/bash"} {if ($10 != $2) {if ($2 != ".") {print "egrep -m 1 \"" $2 "[[:space:]]\" "dbsnp} else {print "echo "}; if ($10 != ".") {print "egrep -m 1 \"" $10 "[[:space:]]\" "dbsnp} else {print "echo "}; print "echo ----------------------"}}' $tmpdir/$vartype.avinput.hg19_snp137_dropped > $tmpdir/$vartype.egreps_dropped.sh
awk -v dbsnp=$dbsnp 'BEGIN{OFS=""; print "#!/bin/bash"} {if ($8 != ".") {print "egrep -m 1 \"" $8 "[[:space:]]\" "dbsnp; print "echo ----------------------"}}' $tmpdir/$vartype.avinput.hg19_snp137_filtered > $tmpdir/$vartype.egreps_filtered.sh
chmod 0700 $tmpdir/$vartype.egreps_dropped.sh $tmpdir/$vartype.egreps_filtered.sh
done
fi

if [ -n "$SNPSIFT" ]
then
for vartype in SNP INDEL
do
java -Xmx1g -Xms1g -jar $SNPSIFT annotate $dbsnp <(gunzip -c $tmpdir/$vartype.vcf.gz|awk 'BEGIN{OFS="\t"} /^#/{print $0} !/^#/{printf("%s\t%s\t.\t",$1,$2); for(i=4;i<=NF;++i)printf("%s\t", $i);printf("\t%s\n", $3)}') 2>$tmpdir/$vartype.snpsift.log > $tmpdir/$vartype.snpsift.vcf &
done
wait
fi

echo "Getting stats on SNPs and indels"
for vartype in SNP INDEL
do
vcf-stats $tmpdir/$vartype".vcf.gz" -p $tmpdir/$vartype &
done

echo "Comparing with NIST calls"
vcf-compare -a $tmpdir/SNP.vcf.gz $NISTVCF &>$tmpdir/SNP.vcf-compare.NIST.txt &
vcf-compare -a -g $tmpdir/SNP.vcf.gz $NISTVCF &>$tmpdir/SNP.vcf-compare.genotyped.NIST.txt &

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
