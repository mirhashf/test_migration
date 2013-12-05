#!/bin/bash

# This script checks for SV deletion events against a truth set and prints out a report for each contig as well as for the whole-genome

set -e

LAKE=/net/kodiak/volumes/lake/shared
RIVER=/net/kodiak/volumes/river/shared

PATH=$PATH:$LAKE/opt/bedtools-2.17.0/bin
REFERENCE=$LAKE/users/marghoob/GATK-bundle-hg19/ucsc.hg19.fa
SYNTHETIC_GENOME=$RIVER/users/marghoob/synthetic_genome/work
SVDELETIONS=$SYNTHETIC_GENOME/deletions.hg19.vcf
TABIX=/usr/lib/bina/tabix/current/bin/tabix

DIR="$( cd "$( dirname "$0" )" && pwd )"

myname=`basename $0`

function usage {
  echo "BREAKDANCER_OUTDIR= CNVNATOR_OUTDIR= BREAKDANCER_OUTDIR= BREAKSEQ_GFF= RECIP_OVERLAP= $myname <workdir> <reportdir>"
  exit 1
}

function print_overlap_counts {
  awk -v mesg="$3" -v min_size=$4 -v max_size=$5 'BEGIN { found = 0; checked = 0; } { size=$3 - $2; if (size < min_size || size >= max_size) next; checked++; if (NF > 4 && $5 != "." && $5 != "0") found++ } END { percent = checked == 0? "nan": 100.0 * found / checked; print found "/" checked " (" percent " %) " mesg }' $1 >> $2
}

WORKDIR=$1
REPORTDIR=$2
[ -z "$BREAKDANCER_OUTDIR" -a -z "$CNVNATOR_OUTDIR" -a -z "$BREAKSEQ_GFF" -a -z "$PINDEL_OUTDIR" ] && usage
[ -z "$WORKDIR" ] && usage

MIN_SIZE=51
MAX_SIZE=10000000

TOOLS=
[ -n "$BREAKDANCER_OUTDIR" ] && TOOLS="breakdancer"
[ -n "$CNVNATOR_OUTDIR" ] && TOOLS="$TOOLS cnvnator"
[ -n "$BREAKSEQ_GFF" ] && TOOLS="$TOOLS breakseq"
[ -n "$PINDEL_OUTDIR" ] && TOOLS="$TOOLS pindel"


mkdir -pv $WORKDIR/truth
for tool in $TOOLS; do
  mkdir -pv $WORKDIR/$tool
  rm -f $WORKDIR/$tool/*
done

if [ -z "$CHR_LIST" ]; then
  CHR_LIST=`awk '{printf $1" "}' $REFERENCE.fai`
fi

REPORT=$WORKDIR/report.txt
rm -f $REPORT $WORKDIR/truth/*

# Convert the SV deletions file to a bedfile
awk -v min_size=$MIN_SIZE -f $DIR/indels.awk <(gunzip -c $SYNTHETIC_GENOME/INDEL.merged.vcf.gz|grep -v SVLEN) | grep deletion | bedtools sort > $WORKDIR/truth/indels.bed
awk -v min_size=$MIN_SIZE -f $DIR/sv_deletions_vcf_to_bed.awk $SVDELETIONS | bedtools sort > $WORKDIR/truth/deletions.bed
cat $WORKDIR/truth/deletions.bed $WORKDIR/truth/indels.bed | bedtools sort | uniq > $WORKDIR/truth/deletions.indels.bed

echo "Generating the BED files"
for chr in $CHR_LIST; do
  echo "Checking results for $chr"
  (
    TRUTH_BED=$WORKDIR/truth/deletions.$chr.bed
    INDEL_BED=$WORKDIR/truth/indels.$chr.bed

    [ -n "$BREAKDANCER_OUTDIR" ] && awk '!/^#/ { if ($7 == "DEL") { print $1"\t"$2"\t"$5"\tBreakdancer" } }' $BREAKDANCER_OUTDIR/$chr.out > $WORKDIR/breakdancer/$chr.bed
    [ -n "$CNVNATOR_OUTDIR" ] && awk '!/^#/ { if ($1 != "deletion") next; split($2, chr_bp_split, ":"); split(chr_bp_split[2], bps, "-"); print chr_bp_split[1]"\t"bps[1]"\t"bps[2]"\tCnvnator" }' $CNVNATOR_OUTDIR/$chr.out > $WORKDIR/cnvnator/$chr.bed
    [ -n "$BREAKSEQ_GFF" ] && (grep "PASS" $BREAKSEQ_GFF | awk -v chr=$chr '{if ($1 == chr && $3 == "Deletion") print $1"\t"$4 - 1 "\t"$5 "\tBreakseq"}' | bedtools sort > $WORKDIR/breakseq/$chr.bed)
    [ -n "$PINDEL_OUTDIR" ] && [ -s "$PINDEL_OUTDIR/$chr.out_D" ] && (grep ChrID $PINDEL_OUTDIR/$chr.out_D | awk -v minsize=$MIN_SIZE '{if ($27 >= 0 && $11 - $10 -1 >= minsize) print $8 "\t" $10 "\t" $11-1 "\tPindel" }' | bedtools sort > $WORKDIR/pindel/$chr.bed)
    awk -v chr=$chr '{if ($1 == chr) print $1"\t"$2 "\t"$3 "\tTruth"}' $WORKDIR/truth/deletions.bed | bedtools sort > $TRUTH_BED
    awk -v chr=$chr '{if ($1 == chr) print $1"\t"$2 "\t"$3 "\tTruth"}' $WORKDIR/truth/deletions.indels.bed | bedtools sort > $INDEL_BED
  ) &
done
wait

echo "Performing comparisons"
for chr in $CHR_LIST; do
  TRUTH_BED=$WORKDIR/truth/deletions.$chr.bed
  INDEL_BED=$WORKDIR/truth/indels.$chr.bed
  for tool in $TOOLS; do
    echo -n "" > $WORKDIR/$tool/$chr.truth.overlap.with.$tool.bed
    echo -n "" > $WORKDIR/$tool/$chr.$tool.overlap.with.truth.bed
  done

  if [ -s "$TRUTH_BED" ]; then
    for tool in $TOOLS; do
      TOOL_BED="$WORKDIR/$tool/$chr.bed"
      if [ -s "$TOOL_BED" ]; then
        bedtools intersect -wao -f $RECIP_OVERLAP -r -a $TRUTH_BED -b $TOOL_BED > $WORKDIR/$tool/$chr.truth.overlap.with.$tool.bed &
      else
        cat $TRUTH_BED > $WORKDIR/$tool/$chr.truth.overlap.with.$tool.bed
      fi
    done
  fi

  if [ -s "$INDEL_BED" ]; then
    for tool in $TOOLS; do
      TOOL_BED="$WORKDIR/$tool/$chr.bed"
      if [ -s "$TOOL_BED" ]; then
        bedtools intersect -wao -f $RECIP_OVERLAP -r -b $INDEL_BED -a $TOOL_BED > $WORKDIR/$tool/$chr.$tool.overlap.with.truth.bed &
      fi
    done
  else
    for tool in $TOOLS; do
      TOOL_BED="$WORKDIR/$tool/$chr.bed"
      [ -s "$TOOL_BED" ] &&  cat $TOOL_BED > $WORKDIR/$tool/$chr.$tool.overlap.with.truth.bed
    done
  fi
done
wait

for tool in $TOOLS; do
  (
    echo "Concatenating stuff for $tool and $CHR_LIST"
    for chr in $CHR_LIST; do
      cat $WORKDIR/$tool/$chr.truth.overlap.with.$tool.bed >> $WORKDIR/$tool/truth.overlap.with.$tool.bed
      cat $WORKDIR/$tool/$chr.$tool.overlap.with.truth.bed >> $WORKDIR/$tool/$tool.overlap.with.truth.bed
    done
    awk -v min_size=$MIN_SIZE '{if (NF > 4 && $5 != "." && $5 != "0" && $3 - $2 >= min_size) print $0}' $WORKDIR/$tool/truth.overlap.with.$tool.bed > $WORKDIR/$tool/truth.found.bed
    awk -v min_size=$MIN_SIZE '{if (NF > 4 && $5 != "." && $5 != "0" && $3 - $2 >= min_size) print $0}' $WORKDIR/$tool/$tool.overlap.with.truth.bed > $WORKDIR/$tool/$tool.found.bed
    awk -v min_size=$MIN_SIZE '{if (!(NF > 4 && $5 != "." && $5 != "0") && $3 - $2 >= min_size) print $1 "\t" $2 "\t" $3 "\t" $4}' $WORKDIR/$tool/truth.overlap.with.$tool.bed > $WORKDIR/$tool/truth.missing.bed
    awk -v min_size=$MIN_SIZE '{if (!(NF > 4 && $5 != "." && $5 != "0") && $3 - $2 >= min_size) print $1 "\t" $2 "\t" $3 "\t" $4}' $WORKDIR/$tool/$tool.overlap.with.truth.bed > $WORKDIR/$tool/$tool.missing.bed
  ) & 
done
wait

echo "Generating the report now"

for chr in $CHR_LIST; do
  for tool in $TOOLS; do
    echo -n ""
    rm -f $WORKDIR/$tool/$chr.*
  done
done

BINS="50,100,200,400,1000,2000,4000,8000,2000000"
for tool in $TOOLS; do
  for prefix in "truth.found" "truth.missing" "$tool.found" "$tool.missing"; do
    print_size_histogram $WORKDIR/$tool/$prefix.bed "$BINS" $WORKDIR/$tool/$prefix
  done 
done

echo "=====================================================================================================" >> $REPORT
echo "all" >> $REPORT

for tool in $TOOLS; do
  print_overlap_counts $WORKDIR/$tool/truth.overlap.with.$tool.bed $REPORT "truth events found in $tool output" $MIN_SIZE $MAX_SIZE
  print_overlap_counts $WORKDIR/$tool/$tool.overlap.with.truth.bed $REPORT "$tool events found in truth set" $MIN_SIZE $MAX_SIZE

  echo "" >> $REPORT
  print_overlap_counts_histogram $WORKDIR/$tool/truth.overlap.with.$tool.bed $REPORT "truth events found in $tool output" 0
  print_overlap_counts_histogram $WORKDIR/$tool/$tool.overlap.with.truth.bed $REPORT "$tool events not found in truth set" 1
  echo "" >> $REPORT
done

echo "#Sensitivity" > $WORKDIR/report.csv
for tool in $TOOLS; do
  print_overlap_counts_histogram $WORKDIR/$tool/truth.overlap.with.$tool.bed $REPORT "truth events found in $tool output" 0 $WORKDIR/report.csv $tool
done

echo "" >> $WORKDIR/report.csv
echo "" >> $REPORT
echo "#False Discovery Rate" >> $WORKDIR/report.csv
for tool in $TOOLS; do
  print_overlap_counts_histogram $WORKDIR/$tool/$tool.overlap.with.truth.bed $REPORT "$tool events not found in truth set" 1 $WORKDIR/report.csv $tool $WORKDIR/counts.csv
done
echo "" >> $WORKDIR/report.csv
