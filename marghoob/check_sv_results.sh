#!/bin/bash

# This script checks for SV deletion events against a truth set and prints out a report for each contig as well as for the whole-genome

set -e

PATH=$PATH:~/lake/opt/bedtools-2.17.0/bin
REFERENCE=~/lake/users/marghoob/GATK-bundle-hg19/ucsc.hg19.fa
SVDELETIONS=work/deletions.hg19.vcf
TABIX=/usr/lib/bina/tabix/current/bin/tabix

DIR="$( cd "$( dirname "$0" )" && pwd )"

myname=`basename $0`

function usage {
  echo "BREAKDANCER_OUTDIR=<breakdancer-dir> CNVNATOR_OUTDIR=<cnvnator-dir> WORKDIR=<work> $myname"
  exit 1
}

function print_overlap_counts {
  awk -v mesg="$3" -v min_size=$4 -v max_size=$5 'BEGIN { found = 0; checked = 0; } { size=$3 - $2; if (size < min_size || size >= max_size) next; checked++; if (NF > 4 && $5 != "." && $5 != "0") found++ } END { percent = checked == 0? "nan": 100.0 * found / checked; print found "/" checked " (" percent " %) " mesg }' $1 >> $2
}

function print_size_histogram {
  bedfile=$1
  bins=$2
  prefix=$3
  awk -v prefix=$prefix -v bins="$bins" \
    'function get_bin(sv_size) {
       for (i=0; i <num_bins; i++) {
         if (sv_size <= sizes[i]) { return i; }
       }
       return -1;
     }
     BEGIN {
       num_bins=split(bins, bins_arr, ",") - 1;
       min_size = bins_arr[1]
       sizes[-1] = min_size
       for (i=0; i < num_bins; i++) {
         sizes[i] = bins_arr[i+2]
         counts[i] = 0
       }
     }
     {
       sv_size = $3 - $2;
       if (sv_size <= min_size) next;
       bin_id = get_bin(sv_size)
       if (bin_id < 0) next;
       counts[bin_id]++
     }
     END {
       counts_file=prefix ".counts"
       printf("%d-%d", min_size, sizes[0]) > counts_file
       for (i=1; i<num_bins;i++) {
         printf(",%d-%d", sizes[i-1], sizes[i]) >> counts_file
       }
       printf(",%d-%d\n", min_size, sizes[num_bins-1]) >> counts_file
       total=counts[0]
       printf("%d", counts[0]) >> counts_file
       for (i=1; i < num_bins; i++) {
         printf(",%d", counts[i]) >> counts_file
         total += counts[i]
       }
       printf(",%d", total) >> counts_file

     }' $bedfile
}

function print_overlap_counts_histogram {
  awk -v mesg="$3" -v complement="$4" -v csv_file="$5" -v counts_file="$7" -v tool_name="$6" \
    'function get_bin(sv_size) {
       for (i=0; i <num_bins; i++) {
         if (sv_size > bins[i-1] && sv_size <= bins[i]) { return i; }
       }
       return -1;
    }
    BEGIN {
      num_bins = 12;
      bins[-1] = 50;
      bins[0] = 100; bins[1] = 200; bins[2] = 400; bins[3] = 600; bins[4] = 800; bins[5] = 1000; bins[6] = 2000; bins[7] = 4000; bins[8] = 6000; bins[9] = 10000; bins[10] = 1000000; bins[11] = 1000000000;
      for (i = 0; i < num_bins; i++) {
        checked_count[i] = 0;
        found_count[i] = 0;
      }
      total_checked = 0;
      total_found = 0;
    }
    {
      sv_size=$3 - $2;
      bin_id = get_bin(sv_size);
      if (bin_id == -1) next;
      checked_count[bin_id]++;
      total_checked++;
      if (NF > 4 && $5 != "." && $5 != "0") {
        found_count[bin_id]++;
        total_found++;
      }
    }
    END {
      if (length(csv_file) > 0) printf("%s", tool_name) >> csv_file
      if (length(counts_file) > 0) printf("%s", tool_name) >> counts_file
      total_percent = (total_checked > 0)? (100.0 * total_found / total_checked): "nan";
      if (complement != 0 && total_percent != "nan") total_percent = 100.0 - total_percent;
      for (i = 0; i < num_bins; i++) {
        percent = "nan"
        if (checked_count[i] != 0) {
          percent = 100.0 * found_count[i] / checked_count[i];
          if (complement != 0) percent = 100 - percent
        }
        printf("%d < size <= %d: %d / %d (%g%%) %s\n", bins[i-1], bins[i], found_count[i], checked_count[i], percent, mesg);
        if (percent == "nan") percent = 0
        if (length(csv_file) > 0) printf(",%g", percent) >> csv_file
        if (length(counts_file) > 0) printf(",%d", checked_count[i]) >> counts_file
      }
      if (total_percent == "nan") total_percent = 0
      if (length(csv_file) > 0) printf(",%g", total_percent) >> csv_file
      if (length(csv_file) > 0) printf("\n") >> csv_file
      if (length(counts_file) > 0) printf(",%d", total_checked) >> counts_file
      if (length(counts_file) > 0) printf("\n") >> counts_file
    }' $1 >> $2
} 

[ -z "$BREAKDANCER_OUTDIR" -a -z "$CNVNATOR_OUTDIR" -a -z "$BREAKSEQ_OUTDIR" -a -z "$PINDEL_OUTDIR" ] && usage
[ -z "$WORKDIR" ] && usage

RECIP_OVERLAP=0.5
MIN_SIZE=51
MAX_SIZE=10000000
PRINT_CHR_STATS=false

TOOLS=
[ -n "$BREAKDANCER_OUTDIR" ] && TOOLS="breakdancer"
[ -n "$CNVNATOR_OUTDIR" ] && TOOLS="$TOOLS cnvnator"
[ -n "$BREAKSEQ_OUTDIR" ] && TOOLS="$TOOLS breakseq"
[ -n "$BREAKSEQ2_OUTDIR" ] && TOOLS="$TOOLS breakseq2"
[ -n "$PINDEL_OUTDIR" ] && TOOLS="$TOOLS pindel"

WORKDIR=$PWD/$WORKDIR
LOGDIR=$WORKDIR/logs

mkdir -pv $WORKDIR/truth
for tool in $TOOLS; do
  mkdir -pv $WORKDIR/$tool
  mkdir -pv $LOGDIR/$tool
  rm -f $WORKDIR/$tool/*
done

if [ -z "$CHR_LIST" ]; then
  CHR_LIST=`awk '{printf $1" "}' $REFERENCE.fai`
fi

REPORT=$WORKDIR/report.txt
rm -f $REPORT $WORKDIR/truth/*

# Convert the SV deletions file to a bedfile
ignore_other=0
awk -v ignore_other=$ignore_other -v min_size=$MIN_SIZE -f $DIR/sv_deletions_vcf_to_bed.awk $SVDELETIONS > $WORKDIR/truth/deletions.bed

echo "Generating the BED files"
for chr in $CHR_LIST; do
  echo "Checking results for $chr"
  (
    TRUTH_BED=$WORKDIR/truth/deletions.$chr.bed

    [ -n "$BREAKDANCER_OUTDIR" ] && awk '!/^#/ { if ($7 == "DEL") { print $1"\t"$2"\t"$5"\tBreakdancer" } }' $BREAKDANCER_OUTDIR/$chr.out > $WORKDIR/breakdancer/$chr.bed
    [ -n "$CNVNATOR_OUTDIR" ] && awk '!/^#/ { if ($1 != "deletion") next; split($2, chr_bp_split, ":"); split(chr_bp_split[2], bps, "-"); print chr_bp_split[1]"\t"bps[1]"\t"bps[2]"\tCnvnator" }' $CNVNATOR_OUTDIR/$chr.out > $WORKDIR/cnvnator/$chr.bed
    [ -n "$BREAKSEQ_OUTDIR" ] && (grep "PASS" $BREAKSEQ_OUTDIR/breakseq_out.gff | awk -v chr=$chr '{if ($1 == chr && $3 == "Deletion") print $1"\t"$4 - 1 "\t"$5 - 1 "\tBreakseq"}' | bedtools sort > $WORKDIR/breakseq/$chr.bed)
    [ -n "$BREAKSEQ2_OUTDIR" ] && (grep "PASS" $BREAKSEQ2_OUTDIR/breakseq_out.gff | awk -v chr=$chr '{if ($1 == chr && $3 == "Deletion") print $1"\t"$4 - 1 "\t"$5 - 1 "\tBreakseq"}' | bedtools sort > $WORKDIR/breakseq2/$chr.bed)
    [ -n "$PINDEL_OUTDIR" ] && [ -s "$PINDEL_OUTDIR/$chr.out_D" ] && (grep ChrID $PINDEL_OUTDIR/$chr.out_D | awk -v minsize=$MIN_SIZE '{if ($27 >= 0 && $11 - $10 -1 >= minsize) print $8 "\t" $10 "\t" $11-1 "\tPindel" }' | bedtools sort > $WORKDIR/pindel/$chr.bed)
    awk -v chr=$chr '{if ($1 == chr) print $1"\t"$2 - 1 "\t"$3 - 1 "\tTruth"}' $WORKDIR/truth/deletions.bed | bedtools sort > $TRUTH_BED
  ) &
done
wait

echo "Performing comparisons"
for chr in $CHR_LIST; do
  TRUTH_BED=$WORKDIR/truth/deletions.$chr.bed
  for tool in $TOOLS; do
    echo -n "" > $WORKDIR/$tool/$chr.truth.overlap.with.$tool.bed
    echo -n "" > $WORKDIR/$tool/$chr.$tool.overlap.with.truth.bed
  done

  if [ -s "$TRUTH_BED" ]; then
    for tool in $TOOLS; do
      TOOL_BED="$WORKDIR/$tool/$chr.bed"
      if [ -s "$TOOL_BED" ]; then
        bedtools intersect -wao -f $RECIP_OVERLAP -r -a $TRUTH_BED -b $TOOL_BED > $WORKDIR/$tool/$chr.truth.overlap.with.$tool.bed &
        bedtools intersect -wao -f $RECIP_OVERLAP -r -b <($TABIX -h work/INDEL.merged.vcf.gz $chr) -a $TOOL_BED > $WORKDIR/$tool/$chr.$tool.overlap.with.truth.bed &
      else
        cat $TRUTH_BED > $WORKDIR/$tool/$chr.truth.overlap.with.$tool.bed
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

if [ "$PRINT_CHR_STATS" == "true" ]; then
  echo "Generating per-chromosome stats"
  for chr in $CHR_LIST; do
    echo "=====================================================================================================" >> $REPORT
    echo "$chr" >> $REPORT

    for tool in $TOOLS; do
      [ -s "$WORKDIR/$tool/$chr.truth.overlap.with.$tool.bed" ] && print_overlap_counts $WORKDIR/$tool/$chr.truth.overlap.with.$tool.bed $REPORT "truth events found in $tool output" $MIN_SIZE $MAX_SIZE
      [ -s "$WORKDIR/$tool/$chr.$tool.overlap.with.truth.bed" ] && print_overlap_counts $WORKDIR/$tool/$chr.$tool.overlap.with.truth.bed $REPORT "$tool events found in truth set" $MIN_SIZE $MAX_SIZE
      echo "" >> $REPORT
      [ -s "$WORKDIR/$tool/$chr.truth.overlap.with.$tool.bed" ] && print_overlap_counts_histogram $WORKDIR/$tool/$chr.truth.overlap.with.$tool.bed $REPORT "truth events found in $tool output" 0
      [ -s "$WORKDIR/$tool/$chr.$tool.overlap.with.truth.bed" ] && print_overlap_counts_histogram $WORKDIR/$tool/$chr.$tool.overlap.with.truth.bed $REPORT "$tool events found in truth set" 1
      echo "" >> $REPORT
    done

    echo "=====================================================================================================" >> $REPORT
    echo "" >> $REPORT
  done
fi

for chr in $CHR_LIST; do
  for tool in $TOOLS; do
    rm -f $WORKDIR/$tool/$chr.*
  done
done

BINS="50,100,200,400,600,800,1000,2000,4000,6000,10000,100000000"
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
