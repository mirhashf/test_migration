#!/bin/bash

myname=`basename $0`

function usage {
  echo "$myname <srcid> <dstid> <inputvcf> <outputvcf>"
  echo "Example usage: $myname b36 b37 variants.b36.vcf variants.b37.vcf"
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


[ $# -ne 4 ] && usage

src=$1
dst=$2
inputvcf=$3
outputvcf=$4

refids=("b36" "b37" "hg18" "hg19")
[ $(contains "${refids[@]}" $src) == "n" ] && echo "Src refname $src non-existent" && exit 1
[ $(contains "${refids[@]}" $dst) == "n" ] && echo "Dst refname $dst non-existent" && exit 1

LAKE=/net/kodiak/volumes/lake/shared
export JAVA_HOME=$LAKE/opt/jdk1.7.0_25/
export PATH=$LAKE/opt/vcftools_0.1.11/bin:$JAVA_HOME/bin:$PATH
export PERL5LIB=$LAKE/opt/vcftools_0.1.11/perl
export GATK_DEV=$LAKE/users/marghoob/gatk-protected/
export CHAIN_DIR=$LAKE/users/marghoob/ftp.broadinstitute.org/Liftover_Chain_Files/
export b36=$LAKE/users/marghoob/ftp.broadinstitute.org/bundle/2.5/b36/human_b36_both
export b37=$LAKE/users/marghoob/ftp.broadinstitute.org/bundle/2.5/b37/human_g1k_v37
export hg19=$LAKE/users/marghoob/GATK-bundle-hg19/ucsc.hg19

CHAIN=$CHAIN_DIR/"$src"to"$dst".chain

[ ! -f "$CHAIN" ] && echo "Chain file $CHAIN doesn't exist" && exit 1

$GATK_DEV/public/perl/liftOverVCF.pl -vcf $inputvcf -chain $CHAIN -out $outputvcf -gatk $GATK_DEV -newRef ${!dst} -oldRef ${!src} -tmp $PWD
