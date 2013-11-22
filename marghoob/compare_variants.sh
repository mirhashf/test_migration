#!/bin/bash

function usage {
  echo "REFCALLS= compare_variants.sh <job1> <job2> <workdir> <summaryfile>"
  exit 1
}

export PERL5LIB=$HOME/vcftools_0.1.11/lib/perl5/site_perl
export JAVA_HOME=$HOME/lake/opt/jdk1.7.0_25/
export PATH=$HOME/vcftools_0.1.11/bin:$JAVA_HOME/bin:$PATH
export GATK_JAR=$HOME/lake/opt/gatk-2.7-2-g6bda569/GenomeAnalysisTK.jar
export SNPSIFT=$HOME/lake/opt/snpEff/SnpSift.jar
dbsnp=$HOME/lake/users/marghoob/GATK-bundle-hg19/dbsnp_137.hg19.vcf
reference=$HOME/lake/users/marghoob/GATK-bundle-hg19/ucsc.hg19.fa

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

set -e

ANNOTATE="true"

[ -z "$1" ] && usage
[ -z "$2" ] && usage
[ -z "$3" ] && usage
[ -z "$4" ] && usage
[ -z "$REFCALLS" ] && REFCALLS=false

dir1=$(print_abs_path $1)
dir2=$(print_abs_path $2)
outdir=$(print_abs_path $3)
outfile=$(print_abs_path $4)

rm -f $outfile

DIR="$( cd "$( dirname "$0" )" && pwd )"

mkdir -pv $outdir

[ -d "$dir1" ] && merge_vcfs $dir1 $reference $outdir/all1.pre_annotated.vcf.gz || (cd $outdir && ln -f -s $dir1 all1.pre_annotated.vcf.gz && ln -f -s $dir1.tbi all1.pre_annotated.vcf.gz.tbi)
[ -d "$dir2" ] && merge_vcfs $dir2 $reference $outdir/all2.pre_annotated.vcf.gz || (cd $outdir && ln -f -s $dir2 all2.pre_annotated.vcf.gz && ln -f -s $dir2.tbi all2.pre_annotated.vcf.gz.tbi)
wait

if [ "$ANNOTATE" == "true" ]
then
echo "Annotating variants"
start_time=$(date +%s)
annotate_vcf $outdir/all1.pre_annotated.vcf.gz $outdir/all1.vcf.gz $outdir/snpsift1.log
annotate_vcf $outdir/all2.pre_annotated.vcf.gz $outdir/all2.vcf.gz $outdir/snpsift2.log
wait
end_time=$(date +%s)
diff_time=$(( $end_time - $start_time ))
echo "Annotation took $diff_time seconds"
else
(mv $outdir/all1.pre_annotated.vcf.gz $outdir/all1.vcf.gz; tabix -f $outdir/all1.vcf.gz)
(mv $outdir/all2.pre_annotated.vcf.gz $outdir/all2.vcf.gz; tabix -f $outdir/all2.vcf.gz)
fi

for filter in ALL PASS; do
  for vartype in SNP INDEL; do
    echo "Separating out $vartype for $filter calls"
    filter_and_select_vcf $outdir/all1.vcf.gz $outdir/$filter.$vartype.1.vcf $vartype $filter $outdir/$filter.$vartype.1.log
    filter_and_select_vcf $outdir/all2.vcf.gz $outdir/$filter.$vartype.2.vcf $vartype $filter $outdir/$filter.$vartype.2.log
  done
done
wait

echo "Comparing the indel and SNP vcfs and getting stats"
for filter in ALL PASS; do
  for vartype in SNP INDEL; do
    vcf-compare $outdir/$filter.$vartype.1.vcf.gz $outdir/$filter.$vartype.2.vcf.gz > $outdir/$filter.$vartype.vcf-compare.txt &
    vcf-stats $outdir/$filter.$vartype.1.vcf.gz -p $outdir/$filter.$vartype.1.stats &
    vcf-stats $outdir/$filter.$vartype.2.vcf.gz -p $outdir/$filter.$vartype.2.stats &
  done
done
wait

echo "Generating known and novel subsets of INDELS and SNPs"

for filter in ALL PASS; do
  for vartype in SNP INDEL; do
    mkdir -pv $outdir/$filter/$vartype

    f1=$outdir/$filter.$vartype.1.vcf.gz
    f2=$outdir/$filter.$vartype.2.vcf.gz

    (vcf-isec $f1 $f2 | bgzip > $outdir/$filter/$vartype/common.vcf.gz; tabix -f $outdir/$filter/$vartype/common.vcf.gz) &
    (vcf-isec -c $f1 $f2 | bgzip > $outdir/$filter/$vartype/1.vcf.gz; tabix -f $outdir/$filter/$vartype/1.vcf.gz) &
    (vcf-isec -c $f2 $f1 | bgzip > $outdir/$filter/$vartype/2.vcf.gz; tabix -f $outdir/$filter/$vartype/2.vcf.gz) &
    wait

    for subset in common 1 2; do
      mkdir -pv $outdir/$filter/$vartype/$subset
      (cat <(gunzip -c $outdir/$filter/$vartype/$subset.vcf.gz|grep "^#") <(gunzip -c $outdir/$filter/$vartype/$subset.vcf.gz |grep -v "^#"|awk '{if ($3 != ".") print $0}') | bgzip > $outdir/$filter/$vartype/$subset/known.vcf.gz; tabix -f $outdir/$filter/$vartype/$subset/known.vcf.gz) &
      (cat <(gunzip -c $outdir/$filter/$vartype/$subset.vcf.gz|grep "^#") <(gunzip -c $outdir/$filter/$vartype/$subset.vcf.gz |grep -v "^#"|awk '{if ($3 == ".") print $0}') | bgzip > $outdir/$filter/$vartype/$subset/novel.vcf.gz; tabix -f $outdir/$filter/$vartype/$subset/novel.vcf.gz) &
      wait

      for subsubset in known novel; do
        vcf-stats $outdir/$filter/$vartype/$subset/$subsubset.vcf.gz -p $outdir/$filter/$vartype/$subset/$subsubset &
      done

      vcf-stats $outdir/$filter/$vartype/$subset.vcf.gz -p $outdir/$filter/$vartype/$subset.stats &
    done
  done
done
wait

echo "Getting the different counts"
$DIR/gen_tables.sh $outdir $dir1 $dir2 > $outfile
