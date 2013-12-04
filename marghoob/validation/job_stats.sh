#!/bin/bash

set -e

myname=`basename $0`
function usage {
  echo "$myname <jobdir> <workdir> <reportdir>"
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


[ $# -ne 3 ] && usage
 
jobdir=$1
workdir=$2
reportdir=$3

LAKE=/net/kodiak/volumes/lake/shared

export PERL5LIB=$LAKE/opt/vcftools/perl
export JAVA_HOME=$LAKE/opt/jdk1.7.0_25/
export PATH=$LAKE/opt/vcftools/bin:$LAKE/opt/tabix:$JAVA_HOME/bin:$PATH
export GATK_JAR=$LAKE/opt/gatk-2.7-2-g6bda569/GenomeAnalysisTK.jar
export SNPSIFT=$LAKE/opt/snpEff/SnpSift.jar
export NISTVCF=$LAKE/users/marghoob/NIST/NISThighConf
dbsnp=$LAKE/users/marghoob/GATK-bundle-hg19/dbsnp_137.hg19.vcf
reference=$LAKE/users/marghoob/GATK-bundle-hg19/ucsc.hg19.fa

DIR="$( cd "$( dirname "$0" )" && pwd )"
mkdir -pv $workdir $reportdir

workdir=$(print_abs_path $workdir)
jobdir=$(print_abs_path $jobdir)
reportdir=$(print_abs_path $reportdir)

# This part does the variant-calling stats
merge_vcfs $jobdir/vcfs $reference $workdir/all.pre_annotated.vcf.gz
wait

echo "Annotating variants" >&2
start_time=$(date +%s)
annotate_vcf $workdir/all.pre_annotated.vcf.gz $workdir/all.vcf.gz $workdir/snpsift.log
wait
end_time=$(date +%s)
diff_time=$(( $end_time - $start_time ))
echo "Annotation took $diff_time seconds" >&2

echo "Concatenating the chromosome VCFs" >&2
vcf=$workdir/all.vcf.gz

echo "Generating subsets of NIST" >&2
mkdir -pv $workdir/NIST
for vartype in SNP INDEL; do
  (cd $workdir/NIST && ln -sf $NISTVCF.$vartype.hg19.annotated.vcf.gz $vartype.vcf.gz && ln -sf $NISTVCF.$vartype.hg19.annotated.vcf.gz.tbi $vartype.vcf.gz.tbi)
done

mkdir -pv $workdir/PASS
for vartype in SNP INDEL; do
  echo "Separating out $vartype for PASS calls" >&2
  filter_and_select_vcf $workdir/all.vcf.gz $workdir/PASS/$vartype.vcf $vartype PASS $workdir/PASS/$vartype.log
done
wait

# Count stuff
echo "Counting things" >&2
mkdir -pv $reportdir/counts
for vartype in SNP INDEL; do
  (vcf-compare $workdir/PASS/$vartype.vcf.gz $workdir/NIST/$vartype.vcf.gz | grep ^VN | cut -f 2- | awk -v jobvcf="$workdir/PASS/$vartype.vcf.gz" 'BEGIN {FS="\t"} {if (NF==2) {if (index($2, jobvcf)==1) print "job "$1; else print "NIST " $1;} else print "both " $1}' > $reportdir/counts/PASS.$vartype.tsv) &
  (gunzip -c $workdir/PASS/$vartype.vcf.gz|awk '/^#/ {print $0} !/^#/ {if ($3 != ".") print $0}'|vcf-tstv|awk '{print $1}' > $reportdir/counts/PASS.$vartype.known.tstv) &
  (gunzip -c $workdir/PASS/$vartype.vcf.gz|awk '/^#/ {print $0} !/^#/ {if ($3 == ".") print $0}'|vcf-tstv|awk '{print $1}' > $reportdir/counts/PASS.$vartype.novel.tstv) &
  (gunzip -c $workdir/PASS/$vartype.vcf.gz | java -jar $SNPSIFT extractFields - ID HET HOM| awk 'BEGIN {FS="\t"} (NR>1){id = ($1 == "")? "Novel": "Known"; het = ($2 != "")? "Het": "Hom"; print id "\t" het}'|sort|uniq -c > $reportdir/counts/PASS.$vartype.hetcounts) &
done
wait

# Now get the variant-calling stats
awk_hethom_str="BEGIN {hetcount=0; total=0} {total += \$1; if (\$3 == \"Het\") hetcount += \$1} END {print hetcount * 1.0 / (total - hetcount)}"

SNP_total=`awk 'BEGIN {sum=0} {sum += $1} END {print sum}' $reportdir/counts/PASS.SNP.hetcounts`
SNP_NIST_sensitivity=`awk 'BEGIN {found=0; nist=0} {if ($1 == "both") found = $2; if ($1 == "NIST") nist = $2} END {print 100.0 * found / (found + nist)}' $reportdir/counts/PASS.SNP.tsv`
SNP_known_frac=`awk -v total=$SNP_total 'BEGIN {sum=0} {sum += $1} END {print sum * 100.0 / total}' <(grep Known $reportdir/counts/PASS.SNP.hetcounts)`
SNP_known_tstv=`cat $reportdir/counts/PASS.SNP.known.tstv`
SNP_novel_tstv=`cat $reportdir/counts/PASS.SNP.novel.tstv`
SNP_known_hethom=`awk "$awk_hethom_str" <(grep Known $reportdir/counts/PASS.SNP.hetcounts)`
SNP_novel_hethom=`awk "$awk_hethom_str" <(grep Novel $reportdir/counts/PASS.SNP.hetcounts)`
SNP_hethom=`awk "$awk_hethom_str" $reportdir/counts/PASS.SNP.hetcounts`

INDEL_total=`awk 'BEGIN {sum=0} {sum += $1} END {print sum}' $reportdir/counts/PASS.INDEL.hetcounts`
INDEL_NIST_sensitivity=`awk 'BEGIN {found=0; nist=0} {if ($1 == "both") found = $2; if ($1 == "NIST") nist = $2} END {print 100.0 * found / (found + nist)}' $reportdir/counts/PASS.INDEL.tsv`
INDEL_known_frac=`awk -v total=$INDEL_total 'BEGIN {sum=0} {sum += $1} END {print sum * 100.0 / total}' <(grep Known $reportdir/counts/PASS.INDEL.hetcounts)`
INDEL_known_hethom=`awk "$awk_hethom_str" <(grep Known $reportdir/counts/PASS.INDEL.hetcounts)`
INDEL_novel_hethom=`awk "$awk_hethom_str" <(grep Novel $reportdir/counts/PASS.INDEL.hetcounts)`
INDEL_hethom=`awk "$awk_hethom_str" $reportdir/counts/PASS.INDEL.hetcounts`

# Get some SV stats
breakseq_events=0
breakdancer_events=0
cnvnator_events=0
pindel_events=0

[ -e "$jobdir/breakseq" ] && breakseq_events=`grep PASS $jobdir/breakseq/breakseq.gff|wc -l`
[ -e "$jobdir/breakdancer" ] && breakdancer_events=`cat $jobdir/breakdancer/*.out|grep -v "^#"|wc -l`
[ -e "$jobdir/cnvnator" ] && cnvnator_events=`cat $jobdir/cnvnator/*.out|grep -v "^#"|wc -l`
[ -e "$jobdir/pindel" ] && pindel_events=`cat $jobdir/pindel/*|grep ChrID|wc -l`

# Do some accuracy QC
is_wes=`grep -c bedfile $jobdir/job.json`

SNP_total_min=3000000
INDEL_total_min=300000
SNP_known_frac_min=95
INDEL_known_frac_min=80
SNP_known_tstv_min=2
SNP_known_hethom_min=1.2
breakseq_min=100
breakdancer_min=10000
cnvnator_min=5000
pindel_min=500000
if [ "$is_wes" != "0" ]; then
  SNP_total_min=20000
  INDEL_total_min=1000
  breakseq_min=0
  breakdancer_min=0
  cnvnator_min=0
  pindel_min=0
fi

SNP_total_pass=$(echo "$SNP_total >= $SNP_total_min" | bc -l)
INDEL_total_pass=$(echo "$INDEL_total >= $INDEL_total_min" | bc -l)
SNP_known_frac_pass=$(echo "$SNP_known_frac >= $SNP_known_frac_min" | bc -l)
INDEL_known_frac_pass=$(echo "$INDEL_known_frac >= $INDEL_known_frac_min" | bc -l)
SNP_known_tstv_pass=$(echo "$SNP_known_tstv >= $SNP_known_tstv_min" | bc -l)
SNP_known_hethom_pass=$(echo "$SNP_known_hethom >= $SNP_known_hethom_min" | bc -l)
breakseq_pass=$(echo "$breakseq_events >= $breakseq_min" | bc -l)
breakdancer_pass=$(echo "$breakdancer_events >= $breakdancer_min" | bc -l)
cnvnator_pass=$(echo "$cnvnator_events >= $cnvnator_min" | bc -l)
pindel_pass=$(echo "$pindel_events >= $pindel_min" | bc -l)

pass_fields="SNP_known_tstv_pass SNP_known_hethom_pass breakseq_pass breakdancer_pass cnvnator_pass pindel_pass"
all_pass=1
for pass_field in $pass_fields; do
  [ "${!pass_field}" == "0" ] && all_pass=0
done

# Print all the stats and pass/fail collected
rm -f $reportdir/stats.csv
fields="SNP_total SNP_known_frac SNP_NIST_sensitivity SNP_hethom SNP_known_tstv SNP_novel_tstv SNP_known_hethom SNP_novel_hethom INDEL_total INDEL_known_frac INDEL_NIST_sensitivity INDEL_hethom INDEL_known_hethom INDEL_novel_hethom breakseq_events breakdancer_events cnvnator_events pindel_events SNP_total_pass INDEL_total_pass SNP_known_frac_pass INDEL_known_frac_pass SNP_known_tstv_pass SNP_known_hethom_pass breakseq_pass breakdancer_pass cnvnator_pass pindel_pass all_pass"
echo -n "jobdir" >> $reportdir/stats.csv
for field in $fields; do
  echo -n ,$field >> $reportdir/stats.csv
done
echo "" >> $reportdir/stats.csv

echo -n "$jobdir" >> $reportdir/stats.csv
for field in $fields; do
  echo -n ,${!field} >> $reportdir/stats.csv
done
echo "" >> $reportdir/stats.csv
