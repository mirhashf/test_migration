#!/bin/bash
# Given a list of job dirs and labels, prints out a report on the individual jobs and pairwise comparisons
# Input is assumed to be stdin

set -e

DIR="$( cd "$( dirname "$0" )" && pwd )"

myname=`basename $0`
function usage {
  echo "$myname <workdir> <reportdir>"
  exit 1
}

[ $# -ne 2 ] && usage
workdir=$1
reportdir=$2

declare -A jobs
while read line; do
  jobdir=$(echo $line | awk '{print $1}')
  label=$(echo $line | awk '{print $2}')
  jobs["$jobdir"]=$label

  $DIR/job_stats.sh $jobdir $workdir/$label $reportdir/$label

  echo "$label"
  cat $reportdir/$label/stats.csv
  echo ""
done < "/proc/${$}/fd/0"

for job1dir in "${!jobs[@]}"; do
  for job2dir in "${!jobs[@]}"; do
    [ "$job1dir" == "$job2dir" ] && continue
    [ "$job1dir" \< "$job2dir" ] && continue

    label1=${jobs[$job1dir]}
    label2=${jobs[$job2dir]}
    LABEL1="$label1" LABEL2="$label2" $DIR/compare_jobs.sh "$job1dir" "$job2dir" $workdir/$label1.vs.$label2 $reportdir/$label1.vs.$label2
    echo "$label1 vs. $label2"
    cat $reportdir/$label1.vs.$label2/stats.csv
    echo ""
  done
done
