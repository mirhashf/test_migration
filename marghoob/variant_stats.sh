#!/bin/bash

set -e

myname=`basename $0`
function usage {
  echo "$myname <jobdir> <workdir> <outputfile> <gatk-version>"
  exit 1
}

function print_abs_path {
  echo $(cd $(dirname $1); pwd)/$(basename $1)
}

function merge_vcfs {
  local vcfdir=$1
  local reference=$2
  local outfile=$3

  local file_list=
  for chr in `awk '{ print $1 }' $reference.fai`; do
    [ -e "$vcfdir/$chr.vcf.gz" ] && file_list="$file_list $vcfdir/$chr.vcf.gz"
  done

  echo "Concatening vcfs from $vcfdir"
  (vcf-concat $file_list | bgzip > $outfile; tabix -f $outfile) &
}

function annotate_vcf {
  local invcf=$1
  local outvcf=$2
  local logfile=$3

  (java  -Xmx1g -Xms1g -jar $SNPSIFT annotate $dbsnp <(gunzip -c $invcf|awk 'BEGIN{OFS="\t"} /^#/{print $0} !/^#/{printf("%s\t%s\t.",$1,$2); for(i=4;i<=NF;++i)printf("\t%s", $i);printf("\n")}') 2>$logfile | java  -Xmx1g -Xms1g -jar $SNPSIFT varType - | bgzip > $outvcf; tabix -f $outvcf) &
}

function filter_and_select_vcf {
  local invcf=$1
  local outvcf=$2
  local vartype=$3
  local filter=$4
  local logfile=$5

  local exclude_filter=""
  [ "$filter" == "PASS" ] && exclude_filter="--excludeFiltered"

  (java -Xmx1g -Xms1g -jar $GATK_JAR -T SelectVariants -U LENIENT_VCF_PROCESSING -selectType $vartype $exclude_filter -env -V $invcf -o $outvcf -R $reference &>$logfile; bgzip -f $outvcf; tabix -f $outvcf.gz) &
}


[ $# -ne 4 ] && usage
 
jobdir=$1
workdir=$2
outfile=$3
gatk_ver=$4

LAKE=/net/kodiak/volumes/lake/shared
RIVER=/net/kodiak/volumes/river/shared

declare -A gatk_paths
gatk_paths["1.6"]="$LAKE/opt/gatk-1.6"
gatk_paths["2.3-9"]="$LAKE/opt/gatk-2.3-9-gd785397"
gatk_paths["2.5-2"]="$LAKE/opt/CancerAnalysisPackage-2013.2-18-g8207e53"
gatk_paths["2.6-5"]="$LAKE/opt/gatk-2.6-5-gba531bd"
gatk_paths["2.7-2"]="$LAKE/opt/gatk-2.7-2-g6bda569"

GATK_PATH=${gatk_paths[$gatk_ver]}

export PERL5LIB=$LAKE/opt/vcftools/perl
export JAVA_HOME=$LAKE/opt/jdk1.7.0_25/
export PATH=$LAKE/opt/vcftools/bin:$LAKE/opt/tabix:$JAVA_HOME/bin:$PATH
export GATK_JAR=$GATK_PATH/GenomeAnalysisTK.jar
export SNPSIFT=$LAKE/opt/snpEff/SnpSift.jar
export NISTVCF=$LAKE/users/marghoob/NIST/NISThighConf
dbsnp=$LAKE/users/marghoob/GATK-bundle-hg19/dbsnp_137.hg19.vcf
reference=$LAKE/users/marghoob/GATK-bundle-hg19/ucsc.hg19.fa

ANNOTATE=true
MERGE=true

rm -f $outfile

DIR="$( cd "$( dirname "$0" )" && pwd )"
mkdir -pv $workdir

workdir=$(print_abs_path $workdir)
jobdir=$(print_abs_path $jobdir)
outfile=$(print_abs_path $outfile)

[ -d "$jobdir" ] && merge_vcfs $jobdir $reference $workdir/all.pre_annotated.vcf.gz || (cd $workdir && ln -sf $jobdir all.pre_annotated.vcf.gz && ln -sf $jobdir.tbi all.pre_annotated.vcf.gz.tbi)
wait

if [ "$ANNOTATE" == "true" ]; then
  echo "Annotating variants"
  start_time=$(date +%s)
  annotate_vcf $workdir/all.pre_annotated.vcf.gz $workdir/all.vcf.gz $workdir/snpsift.log
  wait
  end_time=$(date +%s)
  diff_time=$(( $end_time - $start_time ))
  echo "Annotation took $diff_time seconds"
else
  mv $workdir/all.pre_annotated.vcf.gz $workdir/all.vcf.gz && mv 
fi

echo "Concatenating the chromosome VCFs"
vcf=$workdir/all.vcf.gz

if [ -n "$TRUTH_VCF" ]; then
  echo "Generating subsets of truth VCF"
  mkdir -pv $workdir/truth
  for vartype in SNP INDEL; do
    filter_and_select_vcf $TRUTH_VCF $workdir/truth/$vartype.vcf $vartype ALL $workdir/truth/$vartype.log
  done
fi

if [ -n "$CHECK_NIST" ]; then
  echo "Generating subsets of NIST"
  mkdir -pv $workdir/NIST
  for vartype in SNP INDEL; do
    (cd $workdir/NIST && ln -sf $NISTVCF.$vartype.hg19.annotated.vcf.gz $vartype.vcf.gz && ln -sf $NISTVCF.$vartype.hg19.annotated.vcf.gz.tbi $vartype.vcf.gz.tbi)
  done
fi

for filter in ALL PASS; do
  mkdir -pv $workdir/$filter
  for vartype in SNP INDEL; do
    echo "Separating out $vartype for $filter calls"
    filter_and_select_vcf $workdir/all.vcf.gz $workdir/$filter/$vartype.vcf $vartype $filter $workdir/$filter/$vartype.log
  done
done
wait

for filter in ALL PASS; do
  for vartype in SNP INDEL; do
    ( [ -n "$TRUTH_VCF" ] && vcf-compare $workdir/$filter/$vartype.vcf.gz $workdir/truth/$vartype.vcf.gz > $workdir/$filter/$vartype.truth.compare.txt ) &
    ( [ -n "$CHECK_NIST" ] && vcf-compare $workdir/$filter/$vartype.vcf.gz $workdir/NIST/$vartype.vcf.gz > $workdir/$filter/$vartype.NIST.compare.txt ) &
  done
done
wait

exit

if [ -n "$CHECK_NIST" ]; then

  for vartype in SNP INDEL; do
    mkdir -pv $workdir/NIST/$vartype
  done

  echo "Building subsets from NIST calls"
  for vartype in SNP INDEL; do
    (
      echo "Building subsets of NIST $vartype"
      (cd $workdir/NIST && rm -f $vartype.vcf.gz &&  ln -s $NISTVCF.$vartype.hg19.annotated.vcf.gz $vartype.vcf.gz && tabix -f $vartype.vcf.gz)
      (cat <(gunzip -c $workdir/NIST/$vartype.vcf.gz|grep "^#") <(gunzip -c $workdir/NIST/$vartype.vcf.gz|grep -v "^#"|awk '{if ($3 != ".") print $0}') | bgzip > $workdir/NIST/$vartype/known.vcf.gz; tabix -f $workdir/NIST/$vartype/known.vcf.gz) &
      (cat <(gunzip -c $workdir/NIST/$vartype.vcf.gz|grep "^#") <(gunzip -c $workdir/NIST/$vartype.vcf.gz|grep -v "^#"|awk '{if ($3 == ".") print $0}') | bgzip > $workdir/NIST/$vartype/novel.vcf.gz; tabix -f $workdir/NIST/$vartype/novel.vcf.gz) &
    ) &
  done
fi
wait

for vartype in SNP INDEL; do
  echo "Generating subsets ALL, PASS for $vartype"
  (
    (java -Xmx1g -Xms1g -jar $GATK_JAR -T SelectVariants -U LENIENT_VCF_PROCESSING --excludeNonVariants -selectType $vartype -V $workdir/all.annotated.vcf.gz -R $reference -o $workdir/ALL/$vartype.vcf &>$workdir/ALL/SelectVariants.$vartype.log; bgzip -f $workdir/ALL/$vartype.vcf; tabix -f $workdir/ALL/$vartype.vcf.gz) &
    (java -Xmx1g -Xms1g -jar $GATK_JAR -T SelectVariants -U LENIENT_VCF_PROCESSING --excludeNonVariants -selectType $vartype -V $workdir/all.annotated.vcf.gz -R $reference --excludeFiltered -o $workdir/PASS/$vartype.vcf &>$workdir/PASS/SelectVariants.$vartype.log; bgzip -f $workdir/PASS/$vartype.vcf; tabix -f $workdir/PASS/$vartype.vcf.gz) &
    wait
    (vcf-isec -c $workdir/ALL/$vartype.vcf.gz $workdir/PASS/$vartype.vcf.gz | bgzip > $workdir/NONPASS/$vartype.vcf.gz; tabix -f $workdir/NONPASS/$vartype.vcf.gz)
  ) &
done
wait

for filter in ALL PASS; do #ALL PASS NONPASS; do
  for vartype in SNP INDEL; do
    echo "Generating known and novel subsets of subset $filter of $vartype"
    (cat <(gunzip -c $workdir/$filter/$vartype.vcf.gz|grep "^#") <(gunzip -c $workdir/$filter/$vartype.vcf.gz|grep -v "^#"|awk '{if ($3 != ".") print $0}') | bgzip > $workdir/$filter/$vartype/known.vcf.gz; tabix -f $workdir/$filter/$vartype/known.vcf.gz) &
    (cat <(gunzip -c $workdir/$filter/$vartype.vcf.gz|grep "^#") <(gunzip -c $workdir/$filter/$vartype.vcf.gz|grep -v "^#"|awk '{if ($3 == ".") print $0}') | bgzip > $workdir/$filter/$vartype/novel.vcf.gz; tabix -f $workdir/$filter/$vartype/novel.vcf.gz) &
  done
done
wait

for filter in ALL PASS; do #ALL PASS NONPASS; do
  for vartype in SNP INDEL; do
    echo "Generating stats for subset $filter of $vartype"
    vcf-stats $workdir/$filter/$vartype.vcf.gz -p $workdir/$filter/$vartype.stats/ &
    vcf-stats $workdir/$filter/$vartype/known.vcf.gz -p $workdir/$filter/$vartype/known.stats/ &
    vcf-stats $workdir/$filter/$vartype/novel.vcf.gz -p $workdir/$filter/$vartype/novel.stats/ &
  done
done
wait

if [ -n "$CHECK_NIST" ]; then
  echo "Comparing the various sets against NIST high-confidence calls"
  for filter in ALL PASS; do
    for vartype in SNP INDEL; do
      vcf-compare $workdir/$filter/$vartype.vcf.gz $workdir/NIST/$vartype.vcf.gz > $workdir/$filter/vcf-compare.NIST.$vartype.txt &
      if [ -n "$TRUTH_VCF" ]; then
        vcf-compare $workdir/$filter/$vartype.vcf.gz $workdir/truth/$vartype.vcf.gz > $workdir/$filter/vcf-compare.truth.$vartype.txt &
      fi
      for subset in known novel; do
        vcf-compare $workdir/$filter/$vartype/$subset.vcf.gz $workdir/NIST/$vartype/$subset.vcf.gz > $workdir/$filter/$vartype/vcf-compare.NIST.$subset.txt &
      done
    done
  done
fi
wait

function get_counts() {
  local vcf1=$1
  local vcf2=$2
  local counts_file=$3
  local counts=`grep ^VN $counts_file|cut -f 2-|awk -v vcf1="$vcf1" -v vcf2="$vcf2" -F '\t' '{if (NF == 3) {common=$1}; if (NF==2) {if (index($2, vcf1)) {first=$1} else {second=$1}}} END{total1=common+first; total2=common+second; print common","first","second","(100*common/total2)","(100*first/total1)}'`
  echo "$counts,"
  return 0
}

# Now print out the counts
[ -n "$CHECK_NIST" ] && echo "filter,variant,subset,found_in_nist,missing_in_nist,missing_in_this,nist_sens,nist_fdr,titv,fraction" > $outfile
[ -z "$CHECK_NIST" ] && echo "filter,variant,subset,titv,fraction" > $outfile
for filter in ALL PASS; do #ALL PASS NONPASS; do
  for vartype in SNP INDEL; do
    counts=
    [ -n "$CHECK_NIST" ] && counts=$(get_counts $workdir/$filter/$vartype.vcf.gz $workdir/$filter/$vartype.vcf.gz $workdir/$filter/vcf-compare.NIST.$vartype.txt)
    total_base=`grep -v "^#" $workdir/$filter/$vartype.stats/$vartype.stats.counts|head -n 1|awk '{print $1}'`
    titv=`grep -v "^#" $workdir/$filter/$vartype.stats/$vartype.stats.tstv|head -n 1|awk '{print $3}'`
    echo "$filter,$vartype,all,$total_base,$counts$titv,100" >> $outfile

    for subset in known novel; do
      [ -n "$CHECK_NIST" ] && counts=$(get_counts $workdir/$filter/$vartype/$subset.vcf.gz $workdir/NIST/$vartype/$subset.vcf.gz $workdir/$filter/$vartype/vcf-compare.NIST.$subset.txt)
      total=`grep -v "^#" $workdir/$filter/$vartype/$subset.stats/$subset.stats.counts|head -n 1|awk '{print $1}'`
      titv=`grep -v "^#" $workdir/$filter/$vartype/$subset.stats/$subset.stats.tstv|head -n 1|awk '{print $3}'`
      fraction=`echo "scale=2; 100.0*$total/$total_base"|bc -l`
      echo "$filter,$vartype,$subset,$total,$counts$titv,$fraction" >> $outfile
    done
  done
done

exit
echo "" >> $outfile
echo "===================================================================================" >> $outfile
echo "Stats w.r.t. the truth set"
for filter in ALL PASS NONPASS; do
  for vartype in SNP INDEL; do
    counts=$(get_counts $workdir/$filter/$vartype.vcf.gz $workdir/$filter/$vartype.vcf.gz $workdir/$filter/vcf-compare.NIST.$vartype.txt)
  done
done
#rm -rf $workdir
