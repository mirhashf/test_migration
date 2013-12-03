#!/bin/bash

set -ex

myname=`basename $0`
function usage {
  echo "GENOME= START_E= END_E= NLANES= TOTAL_COVERAGE= myname <basedir>"
  echo "START_E=Error rate on the first base."
  echo "END_E=Error rate on the last base."
  echo "NLANES=Number of lanes. NLANES pairs of FASTQs will be generated."
  echo "TOTAL_COVERAGE=Total coverage across all lanes."
  echo "BEDFILE= optional"
  exit 1
}

LAKE=/net/kodiak/volumes/lake/shared
DIR="$( cd "$( dirname "$0" )" && pwd )"

BASEDIR=$1
[ -z "$BASEDIR" ] && usage
[ -z "$START_E" -o -z "$END_E" -o -z "$NLANES" -o -z "$TOTAL_COVERAGE" ] && usage

mkdir -pv $BASEDIR

DWGSIM=$LAKE/opt/dwgsim/dwgsim

FASTQDIR=$BASEDIR/fastq.$START_E"_"$END_E"_"$NLANES"_"$TOTAL_COVERAGE
mkdir -pv $FASTQDIR
COVERAGE_PER_LANE=`echo "scale = 3; $TOTAL_COVERAGE / 2.0 / $NLANES"|bc -l`
NLANES_MINUS=$((NLANES - 1))

rm -f $FASTQDIR/README
for lane in `seq 0 $NLANES_MINUS`; do
  echo "$DWGSIM -e $START_E,$END_E -E $START_E,$END_E -y 0.02 -r 0 -F 0 -d 330 -s 70 -C $COVERAGE_PER_LANE -1 100 -2 100 -z $lane $GENOME $FASTQDIR/simulated.lane$lane" >> $FASTQDIR/README
  $DWGSIM -e $START_E,$END_E -E $START_E,$END_E -y 0.02 -r 0 -F 0 -d 330 -s 70 -C $COVERAGE_PER_LANE -1 100 -2 100 -z $lane $GENOME $FASTQDIR/simulated.lane$lane
  gzip -v -1 -f $FASTQDIR/simulated.lane$lane.bwa.read1.fastq &>$FASTQDIR/$lane.read1.gzip.log &
  gzip -v -1 -f $FASTQDIR/simulated.lane$lane.bwa.read2.fastq &>$FASTQDIR/$lane.read2.gzip.log &
  wait
done

# generate the jsons for different coverages
BEDFILE_OPT=
[ -n "$BEDFILE" ] && BEDFILE_OPT="--bedfile $BASEDIR/../genome/regions.bed"
$DIR/gen_jsons.py --prefix $FASTQDIR/simulated.lane --nlanes $NLANES --coverage_per_lane $COVERAGE_PER_LANE $BEDFILE_OPT > $FASTQDIR/datasets.json
