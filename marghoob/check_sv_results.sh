#!/bin/bash

# This script checks for SV deletion events against a truth set and prints out a report for each contig as well as for the whole-genome

set -e

PATH=$PATH:~/lake/opt/bedtools-2.17.0/bin
DIR="$( cd "$( dirname "$0" )" && pwd )"

myname=`basename $0`

function usage {
  echo "BREAKDANCER_OUTDIR=<breakdancer-dir> CNVNATOR_OUTDIR=<cnvnator-dir> WORKDIR=<work> $myname"
  exit 1
}

function print_overlap_counts {
  awk -v mesg="$3" -v min_size=$4 -v max_size=$5 'BEGIN { found = 0; checked = 0; } { size=$3 - $2; if (size < min_size || size >= max_size) next; checked++; if (NF > 4 && $5 != ".") found++ } END { percent = checked == 0? "nan": 100.0 * found / checked; print found "/" checked " (" percent " %) " mesg }' $1 >> $2
}

function print_overlap_counts_histogram {
  awk -v mesg="$3" 'function get_bin(sv_size) {
                      for (i=0; i <num_bins; i++) {
                          if (sv_size <= bins[i]) { return i; }
                      }
                      return -1;
                    }
                    BEGIN {
                        num_bins = 8;
                        bins[-1] = 0;
                        bins[0] = 100; bins[1] = 200; bins[2] = 400; bins[3] = 600; bins[4] = 800; bins[5] = 1000; bins[6] = 1000000; bins[7] = 1000000000;
                        for (i = 0; i < num_bins; i++) {
                            checked_count[i] = 0;
                            found_count[i] = 0;
                        }
                    }
                    {
                        sv_size=$3 - $2;
                        bin_id = get_bin(sv_size);
                        if (bin_id == -1) next;
                        checked_count[bin_id]++;
                        if (NF > 4 && $5 != ".") {
                            found_count[bin_id]++;
                        }
                    }
                    END {
                        for (i = 0; i < num_bins; i++) {
                            percent = (checked_count[i] == 0)? "nan": 100.0 * found_count[i] / checked_count[i];
                            printf("%d < size <= %d: %d / %d (%g) %s\n", bins[i-1], bins[i], found_count[i], checked_count[i], percent, mesg);
                        }
                    }' $1 >> $2
} 

[ -z "$BREAKDANCER_OUTDIR" -a -z "$CNVNATOR_OUTDIR" ] && usage
[ -z "$WORKDIR" ] && usage

RECIP_OVERLAP=0.5
MIN_SIZE=0
MAX_SIZE=10000000

TOOLS=
[ -n "$BREAKDANCER_OUTDIR" ] && TOOLS="breakdancer"
[ -n "$CNVNATOR_OUTDIR" ] && TOOLS="$TOOLS cnvnator"

WORKDIR=$PWD/$WORKDIR
LOGDIR=$WORKDIR/logs

mkdir -pv $WORKDIR/truth
for tool in $TOOLS; do
  mkdir -pv $WORKDIR/$tool
  mkdir -pv $LOGDIR/$tool
  rm -f $WORKDIR/$tool/*
done

REFERENCE=~/lake/users/marghoob/GATK-bundle-hg19/ucsc.hg19.fa
SVDELETIONS=work/deletions.hg19.vcf

if [ -z "$CHR_LIST" ]; then
  CHR_LIST=`awk '{print $1}' $REFERENCE.fai`
fi

REPORT=$WORKDIR/report.txt
rm -f $REPORT $WORKDIR/truth/*

# Convert the SV deletions file to a bedfile
ignore_other=0
awk -v ignore_other=$ignore_other -f $DIR/sv_deletions_vcf_to_bed.awk $SVDELETIONS > $WORKDIR/truth/deletions.bed

for chr in $CHR_LIST; do
  echo "Checking results for $chr"
  TRUTH_BED=$WORKDIR/truth/deletions.$chr.sorted.bed

  [ -n "$BREAKDANCER_OUTDIR" ] && awk '!/^#/ { if ($7 == "DEL") { print $1"\t"$2"\t"$5"\tBreakdancer" } }' $BREAKDANCER_OUTDIR/$chr.out > $WORKDIR/breakdancer/$chr.bed
  [ -n "$CNVNATOR_OUTDIR" ] && awk '!/^#/ { if ($1 != "deletion") next; split($2, chr_bp_split, ":"); split(chr_bp_split[2], bps, "-"); print chr_bp_split[1]"\t"bps[1]"\t"bps[2]"\tCnvnator" }' $CNVNATOR_OUTDIR/$chr.out > $WORKDIR/cnvnator/$chr.bed
  awk -v chr=$chr '{if ($1 == chr) print $1"\t"$2"\t"$3"\tTruth"}' $WORKDIR/truth/deletions.bed > $WORKDIR/truth/deletions.$chr.bed
  [ -s "$WORKDIR/truth/deletions.$chr.bed" ] && bedtools sort -i $WORKDIR/truth/deletions.$chr.bed > $TRUTH_BED

  echo "=====================================================================================================" >> $REPORT
  echo "$chr" >> $REPORT

  for tool in $TOOLS; do
    echo -n "" > $WORKDIR/$tool/$chr.truth.overlap.with.$tool.bed
    echo -n "" > $WORKDIR/$tool/$chr.$tool.overlap.with.truth.bed
  done

  if [ -s "$TRUTH_BED" ]; then
    cat $TRUTH_BED >> $WORKDIR/truth/deletions.sorted.bed

    for tool in $TOOLS; do
      TOOL_BED="$WORKDIR/$tool/$chr.bed"
      if [ -s "$TOOL_BED" ]; then
        bedtools intersect -wao -f $RECIP_OVERLAP -r -a $TRUTH_BED -b $TOOL_BED > $WORKDIR/$tool/$chr.truth.overlap.with.$tool.bed
        bedtools intersect -wao -f $RECIP_OVERLAP -r -b $TRUTH_BED -a $TOOL_BED > $WORKDIR/$tool/$chr.$tool.overlap.with.truth.bed
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

  for tool in $TOOLS; do
    cat $WORKDIR/$tool/$chr.truth.overlap.with.$tool.bed >> $WORKDIR/$tool/all.truth.overlap.with.$tool.bed
    cat $WORKDIR/$tool/$chr.$tool.overlap.with.truth.bed >> $WORKDIR/$tool/all.$tool.overlap.with.truth.bed

    if [ -s "$WORKDIR/$tool/$chr.truth.overlap.with.$tool.bed" -o -s "$WORKDIR/$tool/$chr.$tool.overlap.with.truth.bed" ]; then
      print_overlap_counts $WORKDIR/$tool/$chr.truth.overlap.with.$tool.bed $REPORT "truth events found in $tool output" $MIN_SIZE $MAX_SIZE
      print_overlap_counts $WORKDIR/$tool/$chr.$tool.overlap.with.truth.bed $REPORT "$tool events found in truth set" $MIN_SIZE $MAX_SIZE
    fi

    if [ "$PRINT_CHR_HISTOGRAMS" == "true" ]; then
      echo "" >> $REPORT
      [ -s "$WORKDIR/$tool/$chr.truth.overlap.with.$tool.bed" ] && print_overlap_counts_histogram $WORKDIR/$tool/$chr.truth.overlap.with.$tool.bed $REPORT "truth events found in $tool output"
      [ -s "$WORKDIR/$tool/$chr.$tool.overlap.with.truth.bed" ] && print_overlap_counts_histogram $WORKDIR/$tool/$chr.$tool.overlap.with.truth.bed $REPORT "$tool events found in truth set"
      echo "" >> $REPORT
    fi
  done

  echo "=====================================================================================================" >> $REPORT
  echo "" >> $REPORT
done

echo "=====================================================================================================" >> $REPORT
echo "all" >> $REPORT

for tool in $TOOLS; do
  print_overlap_counts $WORKDIR/$tool/all.truth.overlap.with.$tool.bed $REPORT "truth events found in $tool output" $MIN_SIZE $MAX_SIZE
  print_overlap_counts $WORKDIR/$tool/all.$tool.overlap.with.truth.bed $REPORT "$tool events found in truth set" $MIN_SIZE $MAX_SIZE

  print_overlap_counts_histogram $WORKDIR/$tool/all.truth.overlap.with.$tool.bed $REPORT "truth events found in $tool output"
  print_overlap_counts_histogram $WORKDIR/$tool/all.$tool.overlap.with.truth.bed $REPORT "$tool events found in truth set"
done

#echo "=====================================================================================================" >> $REPORT
#echo "Size-based profile of SVs" >> $REPORT
#for tool in $TOOLS; do
#  print_overlap_counts_histogram $WORKDIR/$tool/all.truth.overlap.with.$tool.bed $REPORT "truth events found in $tool output"
#  print_overlap_counts_histogram $WORKDIR/$tool/all.$tool.overlap.with.truth.bed $REPORT "$tool events found in truth set"
#  echo "" >> $REPORT
#done
