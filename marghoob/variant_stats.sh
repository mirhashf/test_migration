#!/bin/bash

set -e

myname=`basename $0`
function usage {
  echo "$myname <jobdir> <workdir> <outputfile> <gatk-version>"
  exit 1
}

[ $# -ne 4 ] && usage
 
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
LAKE_PATH=$HOME/lake

export PERL5LIB=$HOME/vcftools_0.1.11/perl
export JAVA_HOME=$LAKE_PATH/opt/jdk1.7.0_25/
export PATH=$HOME/vcftools_0.1.11/bin:$JAVA_HOME/bin:$PATH
export GATK_JAR=$GATK_PATH/GenomeAnalysisTK.jar
export SNPSIFT=$LAKE_PATH/opt/snpEff/SnpSift.jar
export NISTVCF=$LAKE_PATH/users/marghoob/NIST/NISThighConf
dbsnp=$LAKE_PATH/users/marghoob/GATK-bundle-hg19/dbsnp_137.hg19.vcf
reference=$LAKE_PATH/users/marghoob/GATK-bundle-hg19/ucsc.hg19.fa

ANNOTATE=true

rm -f $outfile

DIR="$( cd "$( dirname "$0" )" && pwd )"
mkdir -pv $tmpdir
awk -v jobdir=$jobdir '{print jobdir"/vcfs/"$1".vcf.gz"}' $reference.fai > $tmpdir/files

echo "Concatenating the chromosome VCFs"
vcf=$tmpdir/all.vcf.gz

if [ "$ANNOTATE" == "true" ]; then
  echo "Clearing ID field and annotating using SnpSift"
  java -Xmx1g -Xms1g -jar $SNPSIFT annotate $dbsnp <(vcf-concat -f $tmpdir/files|awk 'BEGIN{OFS="\t"} /^#/{print $0} !/^#/{printf("%s\t%s\t.",$1,$2); for(i=4;i<=NF;++i)printf("\t%s", $i);printf("\n")}') 2>$tmpdir/snpsift.log | bgzip > $tmpdir/all.annotated.vcf.gz
else
  vcf-concat -f $tmpdir/files | bgzip > $tmpdir/all.annotated.vcf.gz
fi
tabix -f $tmpdir/all.annotated.vcf.gz

for filter in ALL PASS NONPASS; do
  for vartype in SNP INDEL; do
    mkdir -pv $tmpdir/$filter/$vartype
  done
done

for vartype in SNP INDEL; do
  mkdir -pv $tmpdir/NIST/$vartype
done

echo "Building subsets from NIST calls"
for vartype in SNP INDEL; do
  (
    echo "Building subsets of NIST $vartype"
    (cd $tmpdir/NIST && ln -s $NISTVCF.$vartype.hg19.annotated.vcf.gz $vartype.vcf.gz && tabix -f $vartype.vcf.gz)
    (cat <(gunzip -c $tmpdir/NIST/$vartype.vcf.gz|grep "^#") <(gunzip -c $tmpdir/NIST/$vartype.vcf.gz|grep -v "^#"|awk '{if ($3 != ".") print $0}') | bgzip > $tmpdir/NIST/$vartype/known.vcf.gz; tabix -f $tmpdir/NIST/$vartype/known.vcf.gz) &
    (cat <(gunzip -c $tmpdir/NIST/$vartype.vcf.gz|grep "^#") <(gunzip -c $tmpdir/NIST/$vartype.vcf.gz|grep -v "^#"|awk '{if ($3 == ".") print $0}') | bgzip > $tmpdir/NIST/$vartype/novel.vcf.gz; tabix -f $tmpdir/NIST/$vartype/novel.vcf.gz) &
  ) &
done
wait

for vartype in SNP INDEL; do
  echo "Generating subsets ALL, PASS for $vartype"
  (
    (java -Xmx1g -Xms1g -jar $GATK_JAR -T SelectVariants -U LENIENT_VCF_PROCESSING -selectType $vartype -V $tmpdir/all.annotated.vcf.gz -R $reference -o $tmpdir/ALL/$vartype.vcf &>$tmpdir/ALL/SelectVariants.$vartype.log; bgzip -f $tmpdir/ALL/$vartype.vcf; tabix -f $tmpdir/ALL/$vartype.vcf.gz) &
    (java -Xmx1g -Xms1g -jar $GATK_JAR -T SelectVariants -U LENIENT_VCF_PROCESSING -selectType $vartype -V $tmpdir/all.annotated.vcf.gz -R $reference --excludeFiltered -o $tmpdir/PASS/$vartype.vcf &>$tmpdir/PASS/SelectVariants.$vartype.log; bgzip -f $tmpdir/PASS/$vartype.vcf; tabix -f $tmpdir/PASS/$vartype.vcf.gz) &
    wait
    (vcf-isec -c $tmpdir/ALL/$vartype.vcf.gz $tmpdir/PASS/$vartype.vcf.gz | bgzip > $tmpdir/NONPASS/$vartype.vcf.gz; tabix -f $tmpdir/NONPASS/$vartype.vcf.gz)
  ) &
done
wait

for filter in ALL PASS NONPASS; do
  for vartype in SNP INDEL; do
    echo "Generating known and novel subsets of subset $filter of $vartype"
    (cat <(gunzip -c $tmpdir/$filter/$vartype.vcf.gz|grep "^#") <(gunzip -c $tmpdir/$filter/$vartype.vcf.gz|grep -v "^#"|awk '{if ($3 != ".") print $0}') | bgzip > $tmpdir/$filter/$vartype/known.vcf.gz; tabix -f $tmpdir/$filter/$vartype/known.vcf.gz) &
    (cat <(gunzip -c $tmpdir/$filter/$vartype.vcf.gz|grep "^#") <(gunzip -c $tmpdir/$filter/$vartype.vcf.gz|grep -v "^#"|awk '{if ($3 == ".") print $0}') | bgzip > $tmpdir/$filter/$vartype/novel.vcf.gz; tabix -f $tmpdir/$filter/$vartype/novel.vcf.gz) &
  done
done
wait

for filter in ALL PASS NONPASS; do
  for vartype in SNP INDEL; do
    echo "Generating stats for subset $filter of $vartype"
    vcf-stats $tmpdir/$filter/$vartype.vcf.gz -p $tmpdir/$filter/$vartype.stats/ &
    vcf-stats $tmpdir/$filter/$vartype/known.vcf.gz -p $tmpdir/$filter/$vartype/known.stats/ &
    vcf-stats $tmpdir/$filter/$vartype/novel.vcf.gz -p $tmpdir/$filter/$vartype/novel.stats/ &
  done
done
wait

echo "Comparing the various sets against NIST high-confidence calls"
for filter in ALL PASS NONPASS; do
  for vartype in SNP INDEL; do
    vcf-compare $tmpdir/$filter/$vartype.vcf.gz $tmpdir/NIST/$vartype.vcf.gz > $tmpdir/$filter/vcf-compare.NIST.$vartype.txt &
    for subset in known novel; do
      vcf-compare $tmpdir/$filter/$vartype/$subset.vcf.gz $tmpdir/NIST/$vartype/$subset.vcf.gz > $tmpdir/$filter/$vartype/vcf-compare.NIST.$subset.txt &
    done
  done
done
wait

# Now print out the counts
echo "filter,variant,subset,found_in_nist,missing_in_nist,missing_in_this,nist_sens,nist_fdr,titv,fraction" > $outfile
for filter in ALL PASS NONPASS; do
  for vartype in SNP INDEL; do
    vcf1="$tmpdir/$filter/$vartype.vcf.gz"
    vcf2="$tmpdir/NIST/$vartype.vcf.gz"
    counts=`grep ^VN $tmpdir/$filter/vcf-compare.NIST.$vartype.txt|cut -f 2-|awk -v vcf1="$vcf1" -v vcf2="$vcf2" -F '\t' '{if (NF == 3) {common=$1}; if (NF==2) {if (index($2, vcf1)) {first=$1} else {second=$1}}} END{total1=common+first; total2=common+second; print common","first","second","(100*common/total2)","(100*first/total1)}'`
    total_base=`grep -v "^#" $tmpdir/$filter/$vartype.stats/$vartype.stats.counts|head -n 1|awk '{print $1}'`
    titv=`grep -v "^#" $tmpdir/$filter/$vartype.stats/$vartype.stats.tstv|head -n 1|awk '{print $3}'`
    echo "$filter,$vartype,all,$total_base,$counts,$titv,100" >> $outfile

    for subset in known novel; do
      vcf1="$tmpdir/$filter/$vartype/$subset.vcf.gz"
      vcf2="$tmpdir/NIST/$vartype/$subset.vcf.gz"
      counts=`grep ^VN $tmpdir/$filter/$vartype/vcf-compare.NIST.$subset.txt|cut -f 2-|awk -v vcf1="$vcf1" -v vcf2="$vcf2" -F '\t' '{if (NF == 3) {common=$1}; if (NF==2) {if (index($2, vcf1)) {first=$1} else {second=$1}}} END{total1=common+first; total2=common+second; print common","first","second","(100*common/total2)","(100*first/total1)}'`
      total=`grep -v "^#" $tmpdir/$filter/$vartype/$subset.stats/$subset.stats.counts|head -n 1|awk '{print $1}'`
      titv=`grep -v "^#" $tmpdir/$filter/$vartype/$subset.stats/$subset.stats.tstv|head -n 1|awk '{print $3}'`
      fraction=`echo "scale=2; 100.0*$total/$total_base"|bc -l`
      echo "$filter,$vartype,$subset,$total,$counts,$titv,$fraction" >> $outfile
    done
  done
done

#rm -rf $tmpdir
