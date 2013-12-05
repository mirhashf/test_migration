#!/bin/bash

export LAKE=/net/kodiak/volumes/lake/shared

export PERL5LIB=$LAKE/opt/vcftools/perl
export JAVA_HOME=$LAKE/opt/jdk1.7.0_25/
export PATH=$LAKE/opt/vcftools/bin:$LAKE/opt/tabix:$JAVA_HOME/bin:$LAKE/opt/bedtools-2.17.0/bin:$PATH
export GATK_JAR=$LAKE/opt/gatk-2.7-2-g6bda569/GenomeAnalysisTK.jar
export SNPSIFT=$LAKE/opt/snpEff/SnpSift.jar
export NISTVCF=$LAKE/users/marghoob/NIST/NISThighConf
export DBSNP=$LAKE/users/marghoob/GATK-bundle-hg19/dbsnp_137.hg19.vcf
export REFERENCE=$LAKE/users/marghoob/GATK-bundle-hg19/ucsc.hg19.fa

function print_abs_path {
  echo $(cd $(dirname $1); pwd)/$(basename $1)
}

function merge_vcfs {
  local vcfdir=$1
  local reference=$2
  local outfile=$3

  local file_list=
  for chr in `awk '{ print $1 }' $REFERENCE.fai`; do
    [ -e "$vcfdir/$chr.vcf.gz" ] && file_list="$file_list $vcfdir/$chr.vcf.gz"
  done

  echo "Concatening vcfs from $vcfdir" >&2 
  vcf-concat $file_list | bgzip > $outfile
  tabix -f $outfile
}

function annotate_vcf {
  local invcf=$1
  local outvcf=$2
  local logfile=$3

  java  -Xmx1g -Xms1g -jar $SNPSIFT annotate $DBSNP <(gunzip -c $invcf|awk 'BEGIN{OFS="\t"} /^#/{print $0} !/^#/{printf("%s\t%s\t.",$1,$2); for(i=4;i<=NF;++i)printf("\t%s", $i);printf("\n")}') 2>$logfile | java  -Xmx1g -Xms1g -jar $SNPSIFT varType - | bgzip > $outvcf
  tabix -f $outvcf
}

function filter_and_select_vcf {
  local invcf=$1
  local outvcf=$2
  local vartype=$3
  local filter=$4
  local logfile=$5

  local exclude_filter=""
  [ "$filter" == "PASS" ] && exclude_filter="--excludeFiltered"

  java -Xmx1g -Xms1g -jar $GATK_JAR -T SelectVariants -U LENIENT_VCF_PROCESSING -selectType $vartype $exclude_filter -env -V $invcf -o $outvcf -R $REFERENCE &>$logfile
  bgzip -f $outvcf
  tabix -f $outvcf.gz
}

function get_real_basename {
  local path=$1
  
}

export -f print_abs_path
export -f merge_vcfs
export -f annotate_vcf
export -f filter_and_select_vcf
