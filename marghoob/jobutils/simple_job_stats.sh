#!/bin/bash

set -e

myname=`basename $0`
function usage {
  echo "$myname <jobdir> <workdir> <reportdir>"
  exit 1
}

[ $# -ne 3 ] && usage
 
#jobdir=$1
jobvcf=$1
workdir=$2
reportdir=$3

DIR="$( cd "$( dirname "$0" )" && pwd )"

source $DIR/common.sh

mkdir -p $workdir $reportdir

workdir=$(print_abs_path $workdir)
#jobdir=$(print_abs_path $jobdir)
reportdir=$(print_abs_path $reportdir)

#read percent_uniquely_mapped percent_multi_mapped percent_unmapped snp_ti_tv snp_het_hom snp_indel indel_het_hom sequence_quality gc_content duplication_percent <<<"$($DIR/extract_metrics.py $jobdir)"

#read worker_num thread_num time_taken has_avx gatk_version box_type aligner data_type dataset_name <<<"$($DIR/extract_job_json_stats.py $jobdir/job.json)"
echo "$worker_num $thread_num $time_taken $has_avx $gatk_version $box_type $aligner $data_type $dataset_name" >&2

bedfile=
[ "$data_type" == "wes" ] && bedfile=$BEDFILE

if ((1)); then
# This part does the variant-calling stats
echo "Concatenating chromosome VCFs" >&2
#merge_vcfs $jobdir/vcfs $REFERENCE $workdir/all.pre_annotated.vcf.gz
cp $jobvcf $workdir/all.pre_annotated.vcf.gz

echo "Annotating variants" >&2
annotate_vcf $workdir/all.pre_annotated.vcf.gz $workdir/all.vcf.gz $workdir/snpsift.log

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
  (gunzip -c $workdir/PASS/$vartype.vcf.gz|awk '/^#/ {print $0} !/^#/ {if ($3 != ".") print $0}'|vcf-tstv|awk '{print $1}' > $reportdir/counts/PASS.$vartype.known.tstv) &
  (gunzip -c $workdir/PASS/$vartype.vcf.gz|awk '/^#/ {print $0} !/^#/ {if ($3 == ".") print $0}'|vcf-tstv|awk '{print $1}' > $reportdir/counts/PASS.$vartype.novel.tstv) &
  (gunzip -c $workdir/PASS/$vartype.vcf.gz | java -jar $SNPSIFT extractFields - ID HET HOM| awk 'BEGIN {FS="\t"} (NR>1){id = ($1 == "")? "Novel": "Known"; het = ($2 != "")? "Het": "Hom"; print id "\t" het}'|sort|uniq -c > $reportdir/counts/PASS.$vartype.hetcounts) &
done

INS_total=`gunzip -c $workdir/PASS/INDEL.vcf.gz | java -jar $SNPSIFT filter "(exists INS)" | grep -v ^# | grep -v DEL | wc -l`
DEL_total=`gunzip -c $workdir/PASS/INDEL.vcf.gz | java -jar $SNPSIFT filter "(exists DEL)" | grep -v ^# | grep -v INS | wc -l`
wait

INS_DEL_ratio=`echo "scale=2; $INS_total/$DEL_total"|bc -l`

# Now get the variant-calling stats
awk_hethom_str="BEGIN {hetcount=0; total=0} {total += \$1; if (\$3 == \"Het\") hetcount += \$1} END {print hetcount * 1.0 / (total - hetcount)}"

SNP_total=`awk 'BEGIN {sum=0} {sum += $1} END {print sum}' $reportdir/counts/PASS.SNP.hetcounts`
SNP_known_frac=`awk -v total=$SNP_total 'BEGIN {sum=0} {sum += $1} END {print sum * 100.0 / total}' <(grep Known $reportdir/counts/PASS.SNP.hetcounts)`
SNP_known_tstv=`cat $reportdir/counts/PASS.SNP.known.tstv`
SNP_novel_tstv=`cat $reportdir/counts/PASS.SNP.novel.tstv`
SNP_known_hethom=`awk "$awk_hethom_str" <(grep Known $reportdir/counts/PASS.SNP.hetcounts)`
SNP_novel_hethom=`awk "$awk_hethom_str" <(grep Novel $reportdir/counts/PASS.SNP.hetcounts)`
SNP_hethom=`awk "$awk_hethom_str" $reportdir/counts/PASS.SNP.hetcounts`

INDEL_total=`awk 'BEGIN {sum=0} {sum += $1} END {print sum}' $reportdir/counts/PASS.INDEL.hetcounts`
INDEL_known_frac=`awk -v total=$INDEL_total 'BEGIN {sum=0} {sum += $1} END {print sum * 100.0 / total}' <(grep Known $reportdir/counts/PASS.INDEL.hetcounts)`
INDEL_known_hethom=`awk "$awk_hethom_str" <(grep Known $reportdir/counts/PASS.INDEL.hetcounts)`
INDEL_novel_hethom=`awk "$awk_hethom_str" <(grep Novel $reportdir/counts/PASS.INDEL.hetcounts)`
INDEL_hethom=`awk "$awk_hethom_str" $reportdir/counts/PASS.INDEL.hetcounts`

fi

# Get some SV stats
[ -e "$jobdir/breakseq" ] && breakseq_events=`grep PASS $jobdir/breakseq/*.gff|wc -l` || breakseq_events=0
[ -e "$jobdir/breakdancer" ] && breakdancer_events=`cat $jobdir/breakdancer/*.out|grep -v "^#"|wc -l` || breakdancer_events=0
[ -e "$jobdir/cnvnator" ] && cnvnator_events=`cat $jobdir/cnvnator/*.out|grep -v "^#"|wc -l` || cnvnator_events=0
[ -e "$jobdir/pindel" ] && pindel_events=`cat $jobdir/pindel/*|grep ChrID|wc -l` || pindel_events=0

mkdir -pv $workdir/beds
if [ -e "$jobdir/breakdancer" ]; then
  awk '!/^#/ {if ($7 == "DEL" && $8 >= 50) { print $1"\t"$2"\t"$5"\tBreakdancer" }}' $jobdir/breakdancer/*.out > $workdir/beds/breakdancer.bed
  read breakdancer_del breakdancer_del_known_frac breakdancer_sensitivity <<<"$(get_sv_del_stats $workdir/beds/breakdancer.bed)"
fi

if [ -e "$jobdir/breakseq" ]; then
  grep PASS $jobdir/breakseq/*.gff | grep Deletion | awk '{print $1 "\t"$4 - 1 "\t"$5}' | bedtools sort > $workdir/beds/breakseq.bed
  read breakseq_del breakseq_del_known_frac breakseq_sensitivity <<<"$(get_sv_del_stats $workdir/beds/breakseq.bed)"
fi

if [ -e "$jobdir/cnvnator" ]; then
  grep deletion $jobdir/cnvnator/*.out | awk '{split($2, chr_bp_split, ":"); split(chr_bp_split[2], bps, "-"); print chr_bp_split[1]"\t"bps[1]"\t"bps[2]"\tCnvnator"}' | bedtools sort > $workdir/beds/cnvnator.bed
  read cnvnator_del cnvnator_del_known_frac cnvnator_sensitivity <<<"$(get_sv_del_stats $workdir/beds/cnvnator.bed)"
fi

if [ -e "$jobdir/pindel" ]; then
  grep ChrID $jobdir/pindel/*_D | awk '{if ($27 >= 0 && $11 - $10 -1 >= 50) print $8 "\t" $10 "\t" $11-1 "\tPindel" }' | bedtools sort | uniq > $workdir/beds/pindel.bed
  read pindel_del pindel_del_known_frac pindel_sensitivity <<<"$(get_sv_del_stats $workdir/beds/pindel.bed)"
fi

sv_accuracy_fields="breakdancer_del breakdancer_del_known_frac breakseq_del breakseq_del_known_frac cnvnator_del cnvnator_del_known_frac pindel_del pindel_del_known_frac"

# Print all the stats and pass/fail collected
rm -f $reportdir/stats.csv
fields="time_taken  percent_uniquely_mapped percent_multi_mapped percent_unmapped snp_ti_tv snp_het_hom snp_indel indel_het_hom INS_total DEL_total INS_DEL_ratio sequence_quality gc_content duplication_percent SNP_total SNP_known_frac SNP_hethom SNP_known_tstv SNP_novel_tstv SNP_known_hethom SNP_novel_hethom INDEL_total INDEL_known_frac INDEL_hethom INDEL_known_hethom INDEL_novel_hethom breakseq_events breakdancer_events cnvnator_events pindel_events $sv_accuracy_fields"
echo -n "#jobvcf" >> $reportdir/stats.csv
for field in $fields; do
  echo -n ,$field >> $reportdir/stats.csv
done
echo "" >> $reportdir/stats.csv

echo -n "$jobvcf" >> $reportdir/stats.csv
for field in $fields; do
  echo -n ,${!field} >> $reportdir/stats.csv
done
echo "" >> $reportdir/stats.csv
