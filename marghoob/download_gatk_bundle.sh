#!/bin/bash

myname=`basename $0`

function usage {
  echo "$myname <bundleid> <outputdir>"
  echo "<bundleid> must be one of b36,b37,hg18,hg19"
  exit 1
}

# copied from http://stackoverflow.com/questions/3685970/bash-check-if-an-array-contains-a-value
function contains() {
  local n=$#
  local value=${!n}
  for ((i=1;i < $#;i++)) {
    if [ "${!i}" == "${value}" ]; then
      echo "y"
      return 0
    fi
  }
  echo "n"
  return 1
}

[ $# -ne 2 ] && usage

bundleid=$1
outdir=$2

bundleids=("b36" "b37" "hg18" "hg19")

[ $(contains "${bundleids[@]}" $bundleid) == "n" ] && usage

FTP_DIR="ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/2.5/$bundleid"

wget -P $outdir $FTP_DIR/*.fasta.*
wget -P $outdir $FTP_DIR/*.dict.*

for prefix in 1000G dbsnp hapmap Mills
do
wget -P $outdir $FTP_DIR/$prefix*
done

gunzip $outdir/*.gz 
