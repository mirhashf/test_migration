#!/bin/bash

set -e

myname=`basename $0`
function usage {
  echo "$myname <jobdir> <workdir> <reportdir>"
  exit 1
}

[ $# -ne 3 ] && usage
 
jobdir=$1
workdir=$2
reportdir=$3

DIR="$( cd "$( dirname "$0" )" && pwd )"

source $DIR/common.sh

mkdir -p $workdir $reportdir

workdir=$(print_abs_path $workdir)
jobdir=$(print_abs_path $jobdir)
reportdir=$(print_abs_path $reportdir)

# This part does the variant-calling stats
echo "Concatenating chromosome VCFs" >&2
merge_vcfs $jobdir/vcfs $REFERENCE $workdir/all.pre_annotated.vcf.gz

echo "Annotating variants" >&2
annotate_vcf $workdir/all.pre_annotated.vcf.gz $workdir/all.vcf.gz $workdir/snpsift.log

echo "Generating subsets of NIST" >&2
mkdir -p $workdir/NIST
for vartype in SNP INDEL; do
  (cd $workdir/NIST && ln -sf $NISTVCF.$vartype.hg19.annotated.vcf.gz $vartype.vcf.gz && ln -sf $NISTVCF.$vartype.hg19.annotated.vcf.gz.tbi $vartype.vcf.gz.tbi)
done

mkdir -p $workdir/PASS
for vartype in SNP INDEL; do
  echo "Separating out $vartype for PASS calls" >&2
  (filter_and_select_vcf $workdir/all.vcf.gz $workdir/PASS/$vartype.vcf $vartype PASS $workdir/PASS/$vartype.log) &
done
wait

# Count stuff
echo "Counting things" >&2
mkdir -p $reportdir/counts
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
[ -e "$jobdir/breakseq" ] && breakseq_events=`grep PASS $jobdir/breakseq/breakseq.gff|wc -l` || breakseq_events=0
[ -e "$jobdir/breakdancer" ] && breakdancer_events=`cat $jobdir/breakdancer/*.out|grep -v "^#"|wc -l` || breakdancer_events=0
[ -e "$jobdir/cnvnator" ] && cnvnator_events=`cat $jobdir/cnvnator/*.out|grep -v "^#"|wc -l` || cnvnator_events=0
[ -e "$jobdir/pindel" ] && pindel_events=`cat $jobdir/pindel/*|grep ChrID|wc -l` || pindel_events=0

# Do some accuracy QC
is_wes=$(grep -c bedfile $jobdir/job.json|awk '{print $0}')

# Constants ... somewhat arbitrarily chosen
SNP_total_min=3000000
SNP_total_max=5000000
SNP_NIST_sensitivity_min=90
INDEL_total_min=300000
INDEL_total_max=1000000
INDEL_NIST_sensitivity_min=80
SNP_known_frac_min=95
INDEL_known_frac_min=80
SNP_known_tstv_min=2
SNP_known_tstv_max=2.5
SNP_known_hethom_min=1.2
SNP_known_hethom_max=1.8
breakseq_min=0 # aggressive alignment can pretty much make breakseq1 ineffective
breakdancer_min=10000
cnvnator_min=5000
pindel_min=500000
if [ "$is_wes" == "1" ]; then
  SNP_total_min=20000
  SNP_total_max=50000
  SNP_NIST_sensitivity_min=0
  SNP_known_tstv_min=2.5
  SNP_known_tstv_max=3
  INDEL_total_min=1000
  INDEL_total_max=5000
  INDEL_NIST_sensitivity_min=0
  breakseq_min=0
  breakdancer_min=0
  cnvnator_min=0
  pindel_min=0
fi

SNP_total_pass=$(echo "$SNP_total >= $SNP_total_min && $SNP_total <= $SNP_total_max" | bc -l)
INDEL_total_pass=$(echo "$INDEL_total >= $INDEL_total_min && $INDEL_total <= $INDEL_total_max" | bc -l)
SNP_NIST_sensitivity_pass=$(echo "$SNP_NIST_sensitivity >= $SNP_NIST_sensitivity_min" | bc -l)
INDEL_NIST_sensitivity_pass=$(echo "$INDEL_NIST_sensitivity >= $INDEL_NIST_sensitivity_min" | bc -l)
SNP_known_frac_pass=$(echo "$SNP_known_frac >= $SNP_known_frac_min" | bc -l)
INDEL_known_frac_pass=$(echo "$INDEL_known_frac >= $INDEL_known_frac_min" | bc -l)
SNP_known_tstv_pass=$(echo "$SNP_known_tstv >= $SNP_known_tstv_min && $SNP_known_tstv <= $SNP_known_tstv_max" | bc -l)
SNP_known_hethom_pass=$(echo "$SNP_known_hethom >= $SNP_known_hethom_min && $SNP_known_hethom <= $SNP_known_hethom_max" | bc -l)
breakseq_pass=$(echo "$breakseq_events >= $breakseq_min" | bc -l)
breakdancer_pass=$(echo "$breakdancer_events >= $breakdancer_min" | bc -l)
cnvnator_pass=$(echo "$cnvnator_events >= $cnvnator_min" | bc -l)
pindel_pass=$(echo "$pindel_events >= $pindel_min" | bc -l)

pass_fields="SNP_total_pass INDEL_total_pass SNP_NIST_sensitivity_pass INDEL_NIST_sensitivity_pass SNP_known_frac_pass INDEL_known_frac_pass SNP_known_tstv_pass SNP_known_hethom_pass breakseq_pass breakdancer_pass cnvnator_pass pindel_pass"
all_pass=1
for pass_field in $pass_fields; do
  [ "${!pass_field}" == "0" ] && all_pass=0
done

# Print all the stats and pass/fail collected
rm -f $reportdir/stats.csv
fields="SNP_total SNP_known_frac SNP_NIST_sensitivity SNP_hethom SNP_known_tstv SNP_novel_tstv SNP_known_hethom SNP_novel_hethom INDEL_total INDEL_known_frac INDEL_NIST_sensitivity INDEL_hethom INDEL_known_hethom INDEL_novel_hethom breakseq_events breakdancer_events cnvnator_events pindel_events $pass_fields all_pass"
echo -n "#jobdir" >> $reportdir/stats.csv
for field in $fields; do
  echo -n ,$field >> $reportdir/stats.csv
done
echo "" >> $reportdir/stats.csv

echo -n "$jobdir" >> $reportdir/stats.csv
for field in $fields; do
  echo -n ,${!field} >> $reportdir/stats.csv
done
echo "" >> $reportdir/stats.csv
