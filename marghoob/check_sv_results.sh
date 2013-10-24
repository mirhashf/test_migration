#!/bin/bash

set -e

PATH=$PATH:~/lake/opt/bedtools-2.17.0/bin
DIR="$( cd "$( dirname "$0" )" && pwd )"

myname=`basename $0`

function usage {
  echo "BREAKDANCER_OUTDIR=<breakdancer-dir> CNVNATOR_OUTDIR=<cnvnator-dir> WORKDIR=<work> $myname"
  exit 1
}

function print_overlap_counts {
  awk -v mesg="$3" 'BEGIN { found = 0 } { if (NF > 4 && $5 != ".") found++ } END { percent = NR == 0? "nan": 100.0 * found / NR; print found "/" NR " (" percent " %) " mesg }' $1 >> $2
}

[ -z "$BREAKDANCER_OUTDIR" -a -z "$CNVNATOR_OUTDIR" ] && usage
[ -z "$WORKDIR" ] && usage

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
        bedtools intersect -wao -f 0.5 -r -a $TRUTH_BED -b $TOOL_BED > $WORKDIR/$tool/$chr.truth.overlap.with.$tool.bed
        bedtools intersect -wao -f 0.5 -r -b $TRUTH_BED -a $TOOL_BED > $WORKDIR/$tool/$chr.$tool.overlap.with.truth.bed
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

    print_overlap_counts $WORKDIR/$tool/$chr.truth.overlap.with.$tool.bed $REPORT "truth events found in $tool output"
    print_overlap_counts $WORKDIR/$tool/$chr.$tool.overlap.with.truth.bed $REPORT "$tool events found in truth set"
  done

  echo "=====================================================================================================" >> $REPORT
  echo "" >> $REPORT
done

echo "=====================================================================================================" >> $REPORT
echo "all" >> $REPORT

for tool in $TOOLS; do
  print_overlap_counts $WORKDIR/$tool/all.truth.overlap.with.$tool.bed $REPORT "truth events found in $tool output"
  print_overlap_counts $WORKDIR/$tool/all.$tool.overlap.with.truth.bed $REPORT "$tool events found in truth set"
done
