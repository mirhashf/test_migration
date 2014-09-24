#!/bin/bash

export LAKE=/net/kodiak/volumes/lake/shared
export RIVER=/net/kodiak/volumes/river/shared
export DELTA=/net/kodiak/volumes/delta/shared

export PERL5LIB=$LAKE/opt/vcftools/perl
export JAVA_HOME=$LAKE/opt/jdk1.7.0_25/
export PATH=$LAKE/opt/vcftools/bin:$LAKE/opt/tabix:$JAVA_HOME/bin:$LAKE/opt/bedtools-2.17.0/bin:$PATH
export GATK_JAR=$LAKE/opt/GenomeAnalysisTK-2014.2-3.1.7-7-gb0ce0a5/GenomeAnalysisTK.jar
export SNPSIFT=$LAKE/opt/snpEff/SnpSift.jar
#export NISTVCF=$DELTA/ftp/ftp-trace.ncbi.nih.gov/giab/ftp/data/NA12878/variant_calls/NIST/NIST_RTG_PlatGen_merged_highconfidence_v0.2.vcf.gz
export NISTVCF=$DELTA/ftp/ftp-trace.ncbi.nih.gov/giab/ftp/data/NA12878/variant_calls/NIST/NISTIntegratedCalls_13datasets_130719_allcall_UGHapMerge_HetHomVarPASS_VQSRv2.17_all_nouncert_excludesimplerep_excludesegdups_excludedecoy_excludeRepSeqSTRs_noCNVs.vcf.gz
export BEDFILE=$RIVER/users/marghoob/nexterarapidcapture_exome_targetedregions_v1.2.b37.sorted_by_ref.bed
export RTGVCF=$RIVER/users/marghoob/ftp-trace.ncbi.nih.gov/giab/ftp/data/NA12878/variant_calls/RTG/NA12878.phasing_annotated.vcf.gz
export DBSNP=$LAKE/resources/dnaseq/gatk_bundle/2.8/b37/dbsnp_138.b37.vcf
export REFERENCE=$LAKE/references/human/human_g1k_v37_decoy/human_g1k_v37_decoy.fasta
export DGVBED=$RIVER/users/marghoob/git/sandbox/marghoob/svmerge/dgv.bed
export TRUEBED=$RIVER/users/marghoob/svgold/final/gold.merged.bed
#export TRUEBED=$RIVER/users/marghoob/synthetic_genome/work/deletions.b37.bed

export DGVTOTAL=`cat $DGVBED|wc -l`
export TRUETOTAL=`cat $TRUEBED|wc -l`

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

  gunzip -c $invcf | java  -Xmx1g -Xms1g -jar $SNPSIFT varType - | \
    awk '!/^##INFO=<ID=VARTYPE/' | bgzip > $outvcf
  tabix -f $outvcf
}

function filter_and_select_vcf {
  local invcf=$1
  local outvcf=$2
  local vartype=$3
  local filter=$4
  local logfile=$5
  local bedfile=$6

  local exclude_filter=""
  local bedfile_filter=""
  [ "$filter" == "PASS" ] && exclude_filter="--excludeFiltered"
  [ ! -z "$bedfile" ] && bedfile_filter="-L $BEDFILE"

  java -Xmx1g -Xms1g -jar $GATK_JAR -T SelectVariants -U LENIENT_VCF_PROCESSING -selectType $vartype $exclude_filter -env -V $invcf -o $outvcf -R $REFERENCE $bedfile_filter &>$logfile
  bgzip -f $outvcf
  tabix -f $outvcf.gz
}

function get_real_basename {
  local path=$1
  
}

function get_sv_del_stats {
  local sv_bed=$1
  del_count=`cat $sv_bed | wc -l`
  del_count_known=`bedtools intersect -a $sv_bed -b $DGVBED -f 0.5 -u | wc -l`
  del_count_true_found=`bedtools intersect -a $TRUEBED -b $sv_bed -f 0.5 -u | wc -l`
  del_known=0
  del_sensitivity=0
  [ "$del_count" != "0" ] && del_known_frac=`echo "scale=4; 100.0*$del_count_known/$del_count" | bc -l`
  [ "$del_count_true_found" != "0" ] && del_sensitivity=`echo "scale=4; 100.0*$del_count_true_found/$TRUETOTAL" | bc -l`
  echo $del_count $del_known_frac $del_sensitivity
}

export -f print_abs_path
export -f merge_vcfs
export -f annotate_vcf
export -f filter_and_select_vcf
export -f get_sv_del_stats
