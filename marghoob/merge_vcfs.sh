#!/bin/bash

function usage {
  echo "VCFS_DIR=<path to vcfs dir> OUTPUT_VCF=<the merged vcf> BGZIP=<path to bgzip> TABIX=<path to tabix> VCFTOOLS=<path to vcftools dir> REFERENCE=<path to reference> merge_vcfs.sh"
  echo "Make sure the REFERENCE is indexed"
  echo "On a Bina appliance BGZIP=/usr/lib/bina/tabix/current/bin/bgzip TABIX=/usr/lib/bina/tabix/current/bin/tabix"
  echo "vcftools will need to be installed"
  exit 1
}

[ -z "$VCFS_DIR" -o -z "$OUTPUT_VCF" ] && usage

[ -z "$BGZIP" ] && BGZIP=`type -P bgzip`
[ -z "$TABIX" ] && TABIX=`type -P tabix`
[ -z "$VCFTOOLS" ] && VCFCONCAT=`type -P vcf-concat` || VCFCONCAT=$VCFTOOLS/bin/vcf-concat && PERL5LIB=$VCFTOOLS/perl

[ -z "$BGZIP" -o -z "$TABIX" -o -z "$VCFCONCAT" ] && usage

set -ex

file_list=
for chr in `awk '{ print $1 }' $REFERENCE.fai`; do
  [ -e "$VCFS_DIR/$chr.vcf.gz" ] && file_list="$file_list $VCFS_DIR/$chr.vcf.gz"
done

echo "Concatenating vcfs from $VCFS_DIR"
$VCFCONCAT $file_list | $BGZIP > $OUTPUT_VCF.gz 
$TABIX -f $OUTPUT_VCF.gz
