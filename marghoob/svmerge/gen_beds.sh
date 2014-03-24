#!/bin/bash

jobdir=$1
reference="/net/kodiak/volumes/lake/shared/references/human/human_g1k_v37_decoy/human_g1k_v37_decoy.fasta"
reference=b37.fa

function usage {
  echo "script <jobdir>"
  exit 1
}

[ -z "$jobdir" ] && usage

genome=`dirname $jobdir`
study=`dirname $genome`

genome=`basename $genome`
study=`basename $study`
echo $genome $study

beddir=beds/$study/$genome

mkdir -pv $beddir

snpsift="/home/marghoob/lake/opt/snpEff/SnpSift.jar"
gatk="/home/marghoob/lake/opt/CancerAnalysisPackage-2014.1-13-g6b71cb4/GenomeAnalysisTK.jar"
svmerge=$PWD/svmerge.sh

sample_name=`gunzip -c $jobdir/svcnv.vcf.gz | grep ^#CHROM | awk '{print $NF}'`

echo $sample_name

# extract the INDELs from HC
java -Xmx2g -Xms2g -jar $gatk -T SelectVariants -U LENIENT_VCF_PROCESSING -selectType INDEL -V $jobdir/snpindel.reannotated.vcf.gz -o $beddir/indels.vcf -R $reference 
bgzip -c $beddir/indels.vcf > $beddir/indels.marked.vcf.gz

#java -Xmx2g -Xms2g -jar $snpsift vartype $beddir/indels.vcf | bgzip > $beddir/indels.marked.vcf.gz

for vartype in INS DEL; do
  mkdir -pv $beddir/$vartype/raw
  gunzip -c $beddir/indels.marked.vcf.gz | java -Xmx2g -Xms2g -jar $snpsift filter "(exists $vartype)" | \
    grep -v ^# | awk '{size=length($4) - length($5); if ((index($4, ",")==0 && index($5, ",") ==0) && (size>= 50 || size<=-50) && (length($4) == 1 || length($5) == 1)) print $0}' | \
      awk -v vartype=$vartype '{svlen = length($5) - length($4);
                                if (svlen < 0) svlen = -svlen;
                                end = $2 + svlen + 1;
                                if (vartype == "INS") end = $2 + 1;
                                hc_string = "HaplotypeCaller_" $2 "_" end "_" svlen "_0/1"
                                window = (vartype == "INS")? 100: 0;
                                print $1 "\t" $2-window "\t" end+window "\t" hc_string "\t" svlen}' > $beddir/$vartype/raw/HaplotypeCaller.bed
done

gffs=

awk_breakseq_cnvnator="{svlen=(\$5 > 0)? \$5: -\$5; print \$1\"\t\" \$2-1-window \"\t\" \$3+window \"\t\" \$4 \"_\" \$2-1 \"_\" \$3 \"_\" svlen \"_\" \$6 \"\t\" \$5 }"
awk_breakseq="{svlen=(\$5 > 0)? \$5: -\$5; print \$1\"\t\" \$2-1-window \"\t\" \$3+window \"\t\" \$4 \"_\" \$2 \"_\" \$3+1 \"_\" svlen \"_\" \$7 \"_\" \$6 \"\t\" \$5 }"
awk_pindel_breakdancer="{svlen=(\$5 > 0)? \$5: -\$5; print \$1\"\t\" \$2-window \"\t\" \$3+1+window \"\t\" \$4 \"_\" \$2 \"_\" \$3+1 \"_\" svlen \"_\" \$7 \"_\" \$8 \"_\" \$6 \"\t\" \$5 }"
awk_cnvnator="{svlen=(\$5 > 0)? \$5: -\$5; print \$1\"\t\" \$2-1-window \"\t\" \$3+window \"\t\" \$4 \"_\" \$2-1 \"_\" \$3 \"_\" svlen \"_\" \$7 \"_\" \$8 \"_\" \$9 \"_\" \$10 \"_\" \$11 \"_\" \$12 \"_\" \$6 \"\t\" \$5 }"


for svtype in DEL INV INS DUP "DUP:TANDEM"; do
  for tool in Breakdancer BreakSeq Pindel CNVnator; do
    mkdir -pv $beddir/$svtype/raw
    extract_fields="CHROM POS END TOOLNAME SVLEN GEN[0].GT"
    awk_str="$awk_breakseq_cnvnator"
    [ "$tool" == "Pindel" -o "$tool" == "Breakdancer" ] && extract_fields="$extract_fields NORMAL_COUNT NUM_READS"
    [ "$tool" == "CNVnator" ] && extract_fields="$extract_fields natorRD natorP1 natorP2 natorP3 natorP4 natorQ0" && awk_str="$awk_cnvnator"
    [ "$tool" == "Pindel" -o "$tool" == "Breakdancer" ] && awk_str="$awk_pindel_breakdancer"
    window=0
    [ "$tool" == "Pindel" -o "$tool" == "BreakSeq" ] && [ "$svtype" == "INS" ] && window=100
    add_svlen=0
    [ "$tool" == "BreakSeq" -a "$svtype" == "INS" ] && add_svlen=1
    [ "$tool" == "BreakSeq" ] && extract_fields="$extract_fields FILTER" && awk_str="$awk_breakseq"
    gunzip -c $jobdir/svcnv.vcf.gz | \
      java -jar $snpsift filter "(ALT = '<$svtype>') & (TOOLNAME = '$tool')" | \
        awk -v add_svlen=$add_svlen '/^#/ {print $0} !/^#/ {for (i=1; i<=7; i++) printf("%s\t", $i); if (add_svlen==1) printf("SVLEN=100;"); printf("%s", $8); for (i=9; i<=NF; i++) printf("\t%s", $i); printf("\n")}' | \
          java -jar $snpsift filter "(exists SVLEN) & ((SVLEN <= -50) | (SVLEN >= 50))" | \
            java -jar $snpsift extractFields - $extract_fields | grep -v "^#" | \
              awk -v window=$window "$awk_str" | \
                bedtools sort > $beddir/$svtype/raw/$tool.orig.bed
    [ -s $beddir/$svtype/raw/$tool.orig.bed ] && bedtools intersect -a $beddir/$svtype/raw/$tool.orig.bed -b gaps.bed -v > $beddir/$svtype/raw/$tool.bed 
    rm -f $beddir/$svtype/raw/$tool.orig.bed
    [ ! -s $beddir/$svtype/raw/$tool.bed ] && rm -f $beddir/$svtype/raw/$tool.bed
  done
  (cd $beddir/$svtype; $svmerge)
  rm -f $beddir/$svtype/merged2/.merged.f50r.gff
  [ -e "$beddir/$svtype/merged2/.merged.f50r.2.bed" ] && awk -v svtype=$svtype '{print svtype "\t" $0}' $beddir/$svtype/merged2/.merged.f50r.2.bed > $beddir/$svtype/merged2/.merged.f50r.gff
  [ -e "$beddir/$svtype/merged2/.merged.f50r.1.bed" ] && awk -v svtype=$svtype '{print svtype "\t" $0}' $beddir/$svtype/merged2/.merged.f50r.1.bed >> $beddir/$svtype/merged2/.merged.f50r.gff
  gffs="$gffs $beddir/$svtype/merged2/.merged.f50r.gff"
done

set -ex

cat $gffs | ./gen_vcf.py --reference $reference --sample "$sample_name" --sort | bgzip > $beddir/merged.vcf.gz
tabix -f $beddir/merged.vcf.gz

java -jar $gatk -T CombineVariants -R $reference --variant $jobdir/snpindel.vcf.gz --variant $beddir/merged.vcf.gz -o $beddir/snpindel.svcnv.normal.vcf -nt 24 --setKey null
#awk '/^#/ {print} !/^#/ {printf("%s\t%s\t.", $1, $2); for (i=4; i<=NF; i++) {printf("\t%s", $i)}; printf ("\n");}' $beddir/snpindel.svcnv.vcf | java -jar ~/lake/opt/snpEff/SnpSift.jar annotate ~/lake/resources/dnaseq/gatk_bundle/2.8/b37/dbsnp_138.b37.vcf - | bgzip > $beddir/snpindel.svcnv.vcf.gz
bgzip -f $beddir/snpindel.svcnv.normal.vcf
tabix -f $beddir/snpindel.svcnv.normal.vcf.gz

cp $beddir/snpindel.svcnv.normal.vcf.gz* $jobdir
