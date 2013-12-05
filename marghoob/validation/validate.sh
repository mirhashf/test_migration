#!/bin/bash

set -ex

myname=`basename $0`
function usage {
  echo "URL= $myname <basedir>"
  echo "URL example: http://t-rex:19005"
  echo "basedir = Directory to do all stuff in (MUST BE ON RIVER)"
  exit 1
}

BASEDIR=$1
[ -z "$BASEDIR" -o -z "$URL" ] && usage 
DIR="$( cd "$( dirname "$0" )" && pwd )"

mkdir -pv $BASEDIR/output $BASEDIR/jobs

PYTHONPATH="$PYTHONPATH" $DIR/../submission/submit_jobs.py --url "$URL" --datasets $DIR/../real_datasets.json --aligners bwa bina bwamem --enable_vqsr --enable_sv --dataset_names ill_3lane_NA12878 ill_wes_NA12878 --output_dir $BASEDIR/output > $BASEDIR/jobs/jobs.txt

# monitor
cat $BASEDIR/jobs/jobs.txt | $DIR/../submission/monitor.py --url "$URL" > $BASEDIR/jobs/status.txt

# check results
for jobid in `awk '{if ($2 == "Successful") print $1}' $BASEDIR/jobs/status.txt`; do
  $DIR/job_stats.sh $BASEDIR/output/$jobid $BASEDIR/tmp $BASEDIR/reports/$jobid
done
