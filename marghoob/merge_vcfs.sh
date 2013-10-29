#!/bin/bash

function usage {
  echo "VCFS_DIR=<path to vcfs dir> OUTPUT_VCF=<the merged vcf> BGZIP=<path to bgzip> TABIX=<path to tabix> VCFTOOLS=<path to vcftools dir> REFERENCE=<path to reference> merge_vcfs.sh"
  echo "Make sure the REFERENCE is indexed"
  echo "On a Bina appliance BGZIP=/usr/lib/bina/tabix/current/bin/bgzip TABIX=/usr/lib/bina/tabix/current/bin/tabix"
  echo "vcftools will need to be installed"
  exit 1
}

[ -z "$VCFS_DIR" -o -z "$OUTPUT_VCF" -o -z "$BGZIP" -o -z "$VCFTOOLS" -o -z "$VCFTOOLS" ] && usage

set -ex

file_list=
for chr in `awk '{ print $1 }' $REFERENCE.fai`; do
  [ -e "$VCFS_DIR/$chr.vcf.gz" ] && file_list="$file_list $VCFS_DIR/$chr.vcf.gz"
done

echo "Concatenating vcfs from $VCFS_DIR"
export PERL5LIB=$VCFTOOLS/perl
$VCFTOOLS/bin/vcf-concat $file_list | $BGZIP > $OUTPUT_VCF.gz 
$TABIX -f $OUTPUT_VCF.gz
