#!/bin/bash

function usage {
  echo "VCFS_DIR=<path to vcfs dir> OUTPUT_VCF=<the merged vcf> BGZIP=<path to bgzip> TABIX=<path to tabix> REFERENCE=<path to reference> merge_vcfs.sh"
  echo "Make sure the REFERENCE is indexed"
  echo "On a Bina appliance BGZIP=/usr/lib/bina/tabix/current/bin/bgzip TABIX=/usr/lib/bina/tabix/current/bin/tabix"
  echo "vcftools will need to be installed"
  exit 1
}

[ -z "$VCFS_DIR" -o -z "$OUTPUT_VCF" ] && usage

[ -z "$BGZIP" ] && BGZIP=`type -P bgzip`
[ -z "$TABIX" ] && TABIX=`type -P tabix`

[ -z "$BGZIP" -o -z "$TABIX" ] && usage

set -e

file_list=
for chr in `awk '{ print $1 }' $REFERENCE.fai`; do
  vcfgz="$VCFS_DIR/$chr.vcf.gz"
  [ -e "$vcfgz" ] && file_list="$file_list $vcfgz"
done

[ -z "$file_list" ] && echo "No vcfs found" && exit

echo "Concatenating vcfs from $VCFS_DIR"

vcf-concat $file_list | bgzip > $OUTPUT_VCF.gz

$TABIX -f $OUTPUT_VCF.gz
