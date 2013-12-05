#!/bin/bash

set -ex

myname=`basename $0`
function usage {
  echo "CHR_LIST= $myname <truth_vcfs_dir> <jobdir> <report_dir> <tmpdir>"
  exit 1
}

function print_abs_path {
  echo $(cd $(dirname $1); pwd)/$(basename $1)
}

function filter_and_select_vcf {
  local invcf=$1
  local outvcf=$2
  local vartype=$3
  local filter=$4
  local logfile=$5

  local exclude_filter=""
  [ "$filter" == "PASS" ] && exclude_filter="--excludeFiltered"

  (java -Xmx1g -Xms1g -jar $GATK_JAR -T SelectVariants -U LENIENT_VCF_PROCESSING -selectType $vartype $exclude_filter -env -V $invcf -o $outvcf -R $REFERENCE &>$logfile; bgzip -f $outvcf; tabix -f $outvcf.gz) &
}

TRUTH_VCFS_DIR=$1
JOBDIR=$2
REPORTDIR=$3
TMPDIR=$4

[ -z "$TRUTH_VCFS_DIR" -o -z "$JOBDIR" -o -z "$REPORTDIR" -o -z "$TMPDIR" ] && usage

mkdir -pv $REPORTDIR $TMPDIR

LAKE=/net/kodiak/volumes/lake/shared
REFERENCE=$LAKE/users/marghoob/GATK-bundle-hg19/ucsc.hg19.fa
export PERL5LIB=$LAKE/opt/vcftools/perl
export JAVA_HOME=$LAKE/opt/jdk1.7.0_25/
export PATH=$LAKE/opt/vcftools/bin:$LAKE/opt/tabix:$JAVA_HOME/bin:$PATH
export GATK_JAR=$LAKE/opt/gatk-2.7-2-g6bda569/GenomeAnalysisTK.jar
export SNPSIFT=$LAKE/opt/snpEff/SnpSift.jar

mkdir -pv $TMPDIR/truth $TMPDIR/job

TMPDIR=$(print_abs_path $TMPDIR)
TRUTH_VCFS_DIR=$(print_abs_path $TRUTH_VCFS_DIR)

if [ -n "$CHR_LIST" ]; then
  chr_list_to_merge="$CHR_LIST"
  for vartype in SNP INDEL; do
    (rm -f $TMPDIR/truth/$vartype.vcf.gz && tabix -h $TRUTH_VCFS_DIR/$vartype.vcf.gz $chr_list_to_merge | bgzip > $TMPDIR/truth/$vartype.vcf.gz && tabix -f $TMPDIR/truth/$vartype.vcf.gz)
  done
else
  for chromosome in `awk '{print $1}' $REFERENCE.fai`; do
    chr_list_to_merge="$chr_list_to_merge $chromosome"
  done
  for vartype in SNP INDEL; do
    (cd $TMPDIR/truth && ln -sf $TRUTH_VCFS_DIR/$vartype.vcf.gz && ln -sf $TRUTH_VCFS_DIR/$vartype.vcf.gz.tbi)
  done
fi

files_to_merge=
for chromosome in $chr_list_to_merge; do
  [ -e "$JOBDIR/vcfs/$chromosome.vcf.gz" ] && files_to_merge="$files_to_merge $JOBDIR/vcfs/$chromosome.vcf.gz"
done
vcf-concat $files_to_merge | bgzip > $TMPDIR/job/all.vcf.gz; tabix -f $TMPDIR/job/all.vcf.gz

# Now extract the SNPs and INDELs
for vartype in SNP INDEL; do
  filter_and_select_vcf $TMPDIR/job/all.vcf.gz $TMPDIR/job/$vartype.vcf $vartype PASS $TMPDIR/job/SelectVariants.$vartype.log
done
wait

for vartype in SNP INDEL; do
  vcf-compare $TMPDIR/job/$vartype.vcf.gz $TMPDIR/truth/$vartype.vcf.gz > $REPORTDIR/$vartype.vcf-compare.txt
done

# Extract the sensitivity and fdr from the vcf-comparison output
rm -f $REPORTDIR/stats.csv
for vartype in SNP INDEL; do
  grep ^VN $REPORTDIR/$vartype.vcf-compare.txt | cut -f 2- | awk -v vartype=$vartype -v jobvcf="$TMPDIR/job/$vartype.vcf.gz" -v truthvcf="$TMPDIR/truth/$vartype.vcf.gz" 'BEGIN{FS = "\t"; counts[jobvcf]=0; counts[truthvcf] = 0; counts["both"] = 0} {if (NF == 2) {if (index($2, jobvcf) == 1) counts[jobvcf] = $1; else counts[truthvcf] = $1;} else counts["both"] = $1} END {total_true = counts["both"] + counts[truthvcf]; printf("%s,%g,%g\n", vartype, 100.0 * counts["both"] / total_true, 100.0 * counts[jobvcf] / (counts[jobvcf] + counts["both"]))}' >> $REPORTDIR/stats.csv
done

# Now get the stats for the SV tools
[ -e "$JOBDIR/breakdancer" ] && BREAKDANCER_OUTDIR=$JOBDIR/breakdancer
[ -e "$JOBDIR/breakseq" ] && BREAKSEQ_GFF=$JOBDIR/breakseq/breakseq.gff
[ -e "$JOBDIR/cnvnator" ] && CNVNATOR_OUTDIR=$JOBDIR/cnvnator
[ -e "$JOBDIR/pindel" ] && PINDEL_OUTDIR=$JOBDIR/pindel
