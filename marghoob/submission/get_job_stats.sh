#!/bin/bash

set -ex

myname=`basename $0`
function usage {
  echo "CHR_LIST= $myname <jobdir> <reportdir> <tmpdir>"
  exit 1
}

JOBDIR=$
