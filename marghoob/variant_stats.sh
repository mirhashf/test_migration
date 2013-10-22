#!/bin/bash

set -e

myname=`basename $0`
function usage {
  echo "$myname <jobdir> <workdir> <outputfile> <gatk-version>"
  exit 1
}

[ $# -ne 4 ] && usage
 
jobdir=$1
workdir=$2
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
MERGE=true

rm -f $outfile

DIR="$( cd "$( dirname "$0" )" && pwd )"
mkdir -pv $workdir
awk -v jobdir=$jobdir '{print jobdir"/vcfs/"$1".vcf.gz"}' $reference.fai > $workdir/files

echo "Concatenating the chromosome VCFs"
vcf=$workdir/all.vcf.gz

if [ "$ANNOTATE" == "true" ]; then
  echo "Clearing ID field and annotating using SnpSift"
  if [ "$MERGE" == "true" ]; then
    java -Xmx1g -Xms1g -jar $SNPSIFT annotate $dbsnp <(vcf-concat -f $workdir/files|awk 'BEGIN{OFS="\t"} /^#/{print $0} !/^#/{printf("%s\t%s\t.",$1,$2); for(i=4;i<=NF;++i)printf("\t%s", $i);printf("\n")}') 2>$workdir/snpsift.log | bgzip > $workdir/all.annotated.vcf.gz
  else
    java -Xmx1g -Xms1g -jar $SNPSIFT annotate $dbsnp <(gunzip -c $jobdir/all.vcf.gz|awk 'BEGIN{OFS="\t"} /^#/{print $0} !/^#/{printf("%s\t%s\t.",$1,$2); for(i=4;i<=NF;++i)printf("\t%s", $i);printf("\n")}') 2>$workdir/snpsift.log | bgzip > $workdir/all.annotated.vcf.gz
  fi
else
  if [ "$MERGE" == "true" ]; then
    vcf-concat -f $workdir/files | bgzip > $workdir/all.annotated.vcf.gz
  else
    cp $jobdir/all.vcf.gz $workdir/.
  fi
fi
tabix -f $workdir/all.annotated.vcf.gz

for filter in ALL PASS NONPASS; do
  for vartype in SNP INDEL; do
    mkdir -pv $workdir/$filter/$vartype
  done
done

for vartype in SNP INDEL; do
  mkdir -pv $workdir/NIST/$vartype
done

echo "Building subsets from NIST calls"
for vartype in SNP INDEL; do
  (
    echo "Building subsets of NIST $vartype"
    (cd $workdir/NIST && ln -s $NISTVCF.$vartype.hg19.annotated.vcf.gz $vartype.vcf.gz && tabix -f $vartype.vcf.gz)
    (cat <(gunzip -c $workdir/NIST/$vartype.vcf.gz|grep "^#") <(gunzip -c $workdir/NIST/$vartype.vcf.gz|grep -v "^#"|awk '{if ($3 != ".") print $0}') | bgzip > $workdir/NIST/$vartype/known.vcf.gz; tabix -f $workdir/NIST/$vartype/known.vcf.gz) &
    (cat <(gunzip -c $workdir/NIST/$vartype.vcf.gz|grep "^#") <(gunzip -c $workdir/NIST/$vartype.vcf.gz|grep -v "^#"|awk '{if ($3 == ".") print $0}') | bgzip > $workdir/NIST/$vartype/novel.vcf.gz; tabix -f $workdir/NIST/$vartype/novel.vcf.gz) &
  ) &
done
wait

for vartype in SNP INDEL; do
  echo "Generating subsets ALL, PASS for $vartype"
  (
    (java -Xmx1g -Xms1g -jar $GATK_JAR -T SelectVariants -U LENIENT_VCF_PROCESSING -selectType $vartype -V $workdir/all.annotated.vcf.gz -R $reference -o $workdir/ALL/$vartype.vcf &>$workdir/ALL/SelectVariants.$vartype.log; bgzip -f $workdir/ALL/$vartype.vcf; tabix -f $workdir/ALL/$vartype.vcf.gz) &
    (java -Xmx1g -Xms1g -jar $GATK_JAR -T SelectVariants -U LENIENT_VCF_PROCESSING -selectType $vartype -V $workdir/all.annotated.vcf.gz -R $reference --excludeFiltered -o $workdir/PASS/$vartype.vcf &>$workdir/PASS/SelectVariants.$vartype.log; bgzip -f $workdir/PASS/$vartype.vcf; tabix -f $workdir/PASS/$vartype.vcf.gz) &
    wait
    (vcf-isec -c $workdir/ALL/$vartype.vcf.gz $workdir/PASS/$vartype.vcf.gz | bgzip > $workdir/NONPASS/$vartype.vcf.gz; tabix -f $workdir/NONPASS/$vartype.vcf.gz)
  ) &
done
wait

for filter in ALL PASS NONPASS; do
  for vartype in SNP INDEL; do
    echo "Generating known and novel subsets of subset $filter of $vartype"
    (cat <(gunzip -c $workdir/$filter/$vartype.vcf.gz|grep "^#") <(gunzip -c $workdir/$filter/$vartype.vcf.gz|grep -v "^#"|awk '{if ($3 != ".") print $0}') | bgzip > $workdir/$filter/$vartype/known.vcf.gz; tabix -f $workdir/$filter/$vartype/known.vcf.gz) &
    (cat <(gunzip -c $workdir/$filter/$vartype.vcf.gz|grep "^#") <(gunzip -c $workdir/$filter/$vartype.vcf.gz|grep -v "^#"|awk '{if ($3 == ".") print $0}') | bgzip > $workdir/$filter/$vartype/novel.vcf.gz; tabix -f $workdir/$filter/$vartype/novel.vcf.gz) &
  done
done
wait

for filter in ALL PASS NONPASS; do
  for vartype in SNP INDEL; do
    echo "Generating stats for subset $filter of $vartype"
    vcf-stats $workdir/$filter/$vartype.vcf.gz -p $workdir/$filter/$vartype.stats/ &
    vcf-stats $workdir/$filter/$vartype/known.vcf.gz -p $workdir/$filter/$vartype/known.stats/ &
    vcf-stats $workdir/$filter/$vartype/novel.vcf.gz -p $workdir/$filter/$vartype/novel.stats/ &
  done
done
wait

echo "Comparing the various sets against NIST high-confidence calls"
for filter in ALL PASS NONPASS; do
  for vartype in SNP INDEL; do
    vcf-compare $workdir/$filter/$vartype.vcf.gz $workdir/NIST/$vartype.vcf.gz > $workdir/$filter/vcf-compare.NIST.$vartype.txt &
    for subset in known novel; do
      vcf-compare $workdir/$filter/$vartype/$subset.vcf.gz $workdir/NIST/$vartype/$subset.vcf.gz > $workdir/$filter/$vartype/vcf-compare.NIST.$subset.txt &
    done
  done
done
wait

function get_counts() {
  local vcf1=$1
  local vcf2=$2
  local counts_file=$3
  local counts=`grep ^VN $counts_file|cut -f 2-|awk -v vcf1="$vcf1" -v vcf2="$vcf2" -F '\t' '{if (NF == 3) {common=$1}; if (NF==2) {if (index($2, vcf1)) {first=$1} else {second=$1}}} END{total1=common+first; total2=common+second; print common","first","second","(100*common/total2)","(100*first/total1)}'`
  echo $counts
  return 0
}

# Now print out the counts
echo "filter,variant,subset,found_in_nist,missing_in_nist,missing_in_this,nist_sens,nist_fdr,titv,fraction" > $outfile
for filter in ALL PASS NONPASS; do
  for vartype in SNP INDEL; do
    counts=$(get_counts $workdir/$filter/$vartype.vcf.gz $workdir/$filter/$vartype.vcf.gz $workdir/$filter/vcf-compare.NIST.$vartype.txt)
    total_base=`grep -v "^#" $workdir/$filter/$vartype.stats/$vartype.stats.counts|head -n 1|awk '{print $1}'`
    titv=`grep -v "^#" $workdir/$filter/$vartype.stats/$vartype.stats.tstv|head -n 1|awk '{print $3}'`
    echo "$filter,$vartype,all,$total_base,$counts,$titv,100" >> $outfile

    for subset in known novel; do
      counts=$(get_counts $workdir/$filter/$vartype/$subset.vcf.gz $workdir/NIST/$vartype/$subset.vcf.gz $workdir/$filter/$vartype/vcf-compare.NIST.$subset.txt)
      total=`grep -v "^#" $workdir/$filter/$vartype/$subset.stats/$subset.stats.counts|head -n 1|awk '{print $1}'`
      titv=`grep -v "^#" $workdir/$filter/$vartype/$subset.stats/$subset.stats.tstv|head -n 1|awk '{print $3}'`
      fraction=`echo "scale=2; 100.0*$total/$total_base"|bc -l`
      echo "$filter,$vartype,$subset,$total,$counts,$titv,$fraction" >> $outfile
    done
  done
done

#rm -rf $workdir
