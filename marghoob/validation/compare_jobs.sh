#!/bin/bash

set -e

myname=`basename $0`
function usage {
  echo "LABEL1= LABEL2= $myname <job1dir> <job2dir> <workdir> <reportdir>"
  echo "Compare the results of two jobs. Currently only compares variant-calls."
  echo "LABEL1, LABEL2 optional."
  exit 1
}

[ $# -ne 4 ] && usage
 
job1dir=$1
job2dir=$2
workdir=$3
reportdir=$4

DIR="$( cd "$( dirname "$0" )" && pwd )"

source $DIR/common.sh

[ -n "$DEBUG" ] && set -x

mkdir -p $workdir $reportdir

workdir=$(print_abs_path $workdir)
job1dir=$(print_abs_path $job1dir)
job2dir=$(print_abs_path $job2dir)
reportdir=$(print_abs_path $reportdir)

[ -z "$LABEL1" ] && LABEL1=`basename $job1dir`
[ -z "$LABEL2" ] && LABEL2=`basename $job2dir`

if ((0)); then
echo "Concatenating chromosome VCFs" >&2
(merge_vcfs $job1dir/vcfs $REFERENCE $workdir/all1.pre_annotated.vcf.gz) &
(merge_vcfs $job2dir/vcfs $REFERENCE $workdir/all2.pre_annotated.vcf.gz) &
wait

echo "Annotating variants" >&2
(annotate_vcf $workdir/all1.pre_annotated.vcf.gz $workdir/all1.vcf.gz $workdir/snpsift1.log) &
(annotate_vcf $workdir/all2.pre_annotated.vcf.gz $workdir/all2.vcf.gz $workdir/snpsift2.log) &
wait

mkdir -p $workdir/PASS
echo "Separating out INDELs and SNPs" >&2
for vartype in SNP INDEL; do
  (filter_and_select_vcf $workdir/all1.vcf.gz $workdir/PASS/1.$vartype.vcf $vartype PASS $workdir/PASS/1.$vartype.log) &
  (filter_and_select_vcf $workdir/all2.vcf.gz $workdir/PASS/2.$vartype.vcf $vartype PASS $workdir/PASS/2.$vartype.log) &
done
wait

# Generate subsets
echo "Extracting common and private call subsets" >&2
for vartype in SNP INDEL; do
  mkdir -p $workdir/PASS/$vartype
  (vcf-isec $workdir/PASS/1.$vartype.vcf.gz $workdir/PASS/2.$vartype.vcf.gz | bgzip > $workdir/PASS/$vartype/common.vcf.gz; tabix -f $workdir/PASS/$vartype/common.vcf.gz) &
  (vcf-isec -c $workdir/PASS/1.$vartype.vcf.gz $workdir/PASS/2.$vartype.vcf.gz | bgzip > $workdir/PASS/$vartype/1.vcf.gz; tabix -f $workdir/PASS/$vartype/1.vcf.gz) &
  (vcf-isec -c $workdir/PASS/2.$vartype.vcf.gz $workdir/PASS/1.$vartype.vcf.gz | bgzip > $workdir/PASS/$vartype/2.vcf.gz; tabix -f $workdir/PASS/$vartype/2.vcf.gz) &
done
wait

echo "Counting Known, Novel, Het, Hom" >&2
for vartype in SNP INDEL; do
  for subset in 1 2 common; do
    countsdir=$reportdir/counts/PASS/$subset
    mkdir -p $countsdir
    (gunzip -c $workdir/PASS/$vartype/$subset.vcf.gz|awk '/^#/ {print $0} !/^#/ {if ($3 != ".") print $0}'|vcf-tstv|awk '{print $1}' > $countsdir/$vartype.known.tstv) &
    (gunzip -c $workdir/PASS/$vartype/$subset.vcf.gz|awk '/^#/ {print $0} !/^#/ {if ($3 == ".") print $0}'|vcf-tstv|awk '{print $1}' > $countsdir/$vartype.novel.tstv) &
    (gunzip -c $workdir/PASS/$vartype/$subset.vcf.gz | java -jar $SNPSIFT extractFields - ID HET HOM| awk 'BEGIN {FS="\t"} (NR>1){id = ($1 == "")? "Novel": "Known"; het = ($2 != "")? "Het": "Hom"; print id "\t" het}'|sort|uniq -c > $countsdir/$vartype.hetcounts) &
  done
done
fi

echo "Generating BED files from deletion SVs" >&2
for svtool in breakseq breakdancer cnvnator; do
  for i in 1 2; do
    echo -n "" > $workdir/$svtool$i.bed
  done
done

# Do the conversions
[ -e "$job1dir/breakdancer" ] && cat $job1dir/breakdancer/*.out | awk '!/^#/ { if ($7 == "DEL" && $5 - $2 >= 50) { print $1"\t"$2"\t"$5"\tBreakdancer" } }' | bedtools sort > $workdir/breakdancer1.bed
[ -e "$job1dir/cnvnator" ] && cat $job1dir/cnvnator/*.out | awk '!/^#/ { if ($1 != "deletion") next; split($2, chr_bp_split, ":"); split(chr_bp_split[2], bps, "-"); if (bps[2] - bps[1] >= 50) print chr_bp_split[1]"\t"bps[1]"\t"bps[2]"\tCnvnator" }' | bedtools sort > $workdir/cnvnator1.bed
[ -s "$job1dir/breakseq/breakseq.gff" ] && (grep "PASS" $job1dir/breakseq/breakseq.gff | awk '{if ($3 == "Deletion" && $5 - $4 + 1 >= 50) print $1"\t"$4 - 1 "\t"$5 "\tBreakseq"}' | bedtools sort > $workdir/breakseq1.bed)
[ -e "$job1dir/pindel" ] && cat $job1dir/pindel/*._D | grep ChrID | awk -v minsize=50 '{if ($27 >= 0 && $11 - $10 -1 >= minsize) print $8 "\t" $10 "\t" $11-1 "\tPindel" }' | bedtools sort > $workdir/pindel1.bed

[ -e "$job2dir/breakdancer" ] && cat $job2dir/breakdancer/*.out | awk '!/^#/ { if ($7 == "DEL" && $5 - $2 >= 50) { print $1"\t"$2"\t"$5"\tBreakdancer" } }' | bedtools sort > $workdir/breakdancer2.bed
[ -e "$job2dir/cnvnator" ] && cat $job2dir/cnvnator/*.out | awk '!/^#/ { if ($1 != "deletion") next; split($2, chr_bp_split, ":"); split(chr_bp_split[2], bps, "-"); if (bps[2] - bps[1] >= 50) print chr_bp_split[1]"\t"bps[1]"\t"bps[2]"\tCnvnator" }' | bedtools sort > $workdir/cnvnator2.bed
[ -s "$job2dir/breakseq/breakseq.gff" ] && grep "PASS" $job2dir/breakseq/breakseq.gff | awk '{if ($3 == "Deletion" && $5 - $4 + 1 >= 50) print $1"\t"$4 - 1 "\t"$5 "\tBreakseq"}' | bedtools sort > $workdir/breakseq2.bed
[ -e "$job2dir/pindel" ] && cat $job2dir/pindel/*._D | grep ChrID | awk -v minsize=50 '{if ($27 >= 0 && $11 - $10 -1 >= minsize) print $8 "\t" $10 "\t" $11-1 "\tPindel" }' | bedtools sort > $workdir/pindel2.bed

echo "Generating comman and private deletion SV subset counts" >&2
for svtool in breakdancer breakseq cnvnator pindel; do
  mkdir -p $reportdir/counts/PASS/1 $reportdir/counts/PASS/2 $reportdir/counts/PASS/common
  bed1=$workdir/$svtool"1.bed"
  bed2=$workdir/$svtool"2.bed"

  for subset in 1 2 common; do
    echo 0 > $reportdir/counts/PASS/$subset/$svtool.count
  done

  if [ -s "$bed1" ]; then
    if [ -s "$bed2" ]; then
      bedtools intersect -a $bed1 -b $bed2 -r -f 0.5 -u | wc -l > $reportdir/counts/PASS/common/$svtool.count
      bedtools intersect -a $bed1 -b $bed2 -r -f 0.5 -v | wc -l > $reportdir/counts/PASS/1/$svtool.count
      bedtools intersect -b $bed1 -a $bed2 -r -f 0.5 -v | wc -l > $reportdir/counts/PASS/2/$svtool.count
    else
      wc -l $bed1 > $reportdir/counts/PASS/1/$svtool.count
    fi
  else
    if [ -s "$bed2" ]; then
      wc -l $bed2 > $reportdir/counts/PASS/2/$svtool.count
    fi
  fi
done

rm -f $reportdir/stats.csv
awk_hethom_str="BEGIN {hetcount=0; total=0} {total += \$1; if (\$3 == \"Het\") hetcount += \$1} END {if (total - hetcount > 0) print hetcount * 1.0 / (total - hetcount); else print 0}"
awk_total_str="BEGIN {sum=0} {sum += \$1} END {print sum}"

fields="SNP_total INDEL_total SNP_known_frac INDEL_known_frac SNP_hethom INDEL_hethom SNP_known_tstv SNP_novel_tstv SNP_known_hethom INDEL_known_hethom SNP_novel_hethom INDEL_novel_hethom breakseq_count breakdancer_count cnvnator_count pindel_count"
echo -n "#subset" >> $reportdir/stats.csv
for field in $fields; do
  echo -n ,$field >> $reportdir/stats.csv
done
echo "" >> $reportdir/stats.csv

declare -A labels
labels["1"]="$LABEL1-only"
labels["2"]="$LABEL2-only"
labels["common"]="$LABEL1-and-$LABEL2"
for subset in 1 2 common; do
  countsdir=$reportdir/counts/PASS/$subset
  # Now get the variant-calling stats
  SNP_total=`awk "$awk_total_str" $countsdir/SNP.hetcounts`
  SNP_known_frac=`awk -v total=$SNP_total 'BEGIN {sum=0} {sum += $1} END {print sum * 100.0 / total}' <(grep Known $countsdir/SNP.hetcounts)`
  SNP_known_tstv=`cat $countsdir/SNP.known.tstv`
  SNP_novel_tstv=`cat $countsdir/SNP.novel.tstv`
  SNP_known_hethom=`awk "$awk_hethom_str" <(grep Known $countsdir/SNP.hetcounts)`
  SNP_novel_hethom=`awk "$awk_hethom_str" <(grep Novel $countsdir/SNP.hetcounts)`
  SNP_hethom=`awk "$awk_hethom_str" $countsdir/SNP.hetcounts`

  INDEL_total=`awk "$awk_total_str" $countsdir/INDEL.hetcounts`
  INDEL_known_frac=`awk -v total=$INDEL_total 'BEGIN {sum=0} {sum += $1} END {print sum * 100.0 / total}' <(grep Known $countsdir/INDEL.hetcounts)`
  INDEL_known_hethom=`awk "$awk_hethom_str" <(grep Known $countsdir/INDEL.hetcounts)`
  INDEL_novel_hethom=`awk "$awk_hethom_str" <(grep Novel $countsdir/INDEL.hetcounts)`
  INDEL_hethom=`awk "$awk_hethom_str" $countsdir/INDEL.hetcounts`

  breakseq_count=`cat $countsdir/breakseq.count`
  breakdancer_count=`cat $countsdir/breakdancer.count`
  cnvnator_count=`cat $countsdir/cnvnator.count`
  pindel_count=`cat $countsdir/pindel.count`

  # Print all the stats
  echo -n "${labels[$subset]}" >> $reportdir/stats.csv
  for field in $fields; do
    echo -n ,${!field} >> $reportdir/stats.csv
  done
  echo "" >> $reportdir/stats.csv
done
