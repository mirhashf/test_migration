#!/bin/bash

set -ex

myname=`basename $0`
function usage {
  echo "CHR_LIST= START_E= END_E= NLANES= TOTAL_COVERAGE= URL= $myname <basedir>"
  exit 1
}

BASEDIR=$1
[ -z "$BASEDIR" -o -z "$START_E" -o -z "$END_E" -o -z "$NLANES" -o -z "$TOTAL_COVERAGE" -o -z "$URL" ] && usage 
DIR="$( cd "$( dirname "$0" )" && pwd )"

BEDFILE=
[ -n "$CHR_LIST" ] && BEDFILE=$BASEDIR/genome/regions.bed

[ "$SKIP_GENOME_GEN" != "true" ] && CHR_LIST="$CHR_LIST" $DIR/gen_diploid.sh $BASEDIR
[ "$SKIP_FASTQ_GEN" != "true" ] && BEDFILE=$BEDFILE GENOME=$BASEDIR/genome/NA12878.fa START_E=$START_E END_E=$END_E NLANES=$NLANES TOTAL_COVERAGE=$TOTAL_COVERAGE $DIR/gen_reads.sh $BASEDIR/fastqs

# submit
mkdir -pv $BASEDIR/output $BASEDIR/jobs
FASTQDIR=$BASEDIR/fastqs/fastq.$START_E"_"$END_E"_"$NLANES"_"$TOTAL_COVERAGE

if [ "$SKIP_SUBMIT_JOBS" != "true" ]; then
  PYTHONPATH="$PYTHONPATH" $DIR/../submission/submit_jobs.py --url "$URL" --datasets $FASTQDIR/datasets.json --output_dir $BASEDIR/output > $BASEDIR/jobs/jobs.txt

  # monitor
  cat $BASEDIR/jobs/jobs.txt | $DIR/../submission/monitor.py --url "$URL" > $BASEDIR/jobs/status.txt

  # check results
  for jobid in `awk '{if ($2 == "Successful") print $1}' $BASEDIR/jobs/status.txt`; do
    CHR_LIST="$CHR_LIST" $DIR/get_stats.sh $BASEDIR/vcfs $BASEDIR/output/$jobid $BASEDIR/reports/$jobid $BASEDIR/tmp
  done
fi
