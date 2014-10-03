#!/bin/bash

R=$(cd $(dirname $0) && pwd -P)

set -x

export BINA_FRONTEND=t-rex:19011
export BINA_USER=$USER 
export BINA_PASS=b 
FASTA=lake:/references/human/hg19.major/hg19.major.fa
DBSNP=lake:/resources/dnaseq/gatk_bundle/2.8/hg19/dbsnp_138.hg19.vcf
PON=river:/users/marghoob/cosmic/refseq_exome_10bp_hg19_300_1kg_normal_panel.hg19.vcf
COSMIC=river:/users/marghoob/cosmic/cosmic.hg19.vcf
OUTPUT_PRE_PREFIX=river:/users/amir/gel
BAM_PREFIX=delta:/prj/GELA20140221/cancer/raw

run_mutect() {
	$R/submit-mutect-on-bams fasta "$FASTA" dbsnp "$DBSNP" pon "$PON" cosmic "$COSMIC" study "$1" output "$OUTPUT_PRE_PREFIX/$1" tumor "$BAM_PREFIX/$1/$2/${2}_SEQA/$2.bam" normal "$BAM_PREFIX/$1/$3/${3}_SEQA/$3.bam" < $R/nothing.json
}

run_mutect GELS00000000008 GELG00000000044 GELG00000000045
run_mutect GELS00000000007 GELG00000000072 GELG00000000073
run_mutect GELS00000000006 GELG00000000040 GELG00000000041
run_mutect GELS00000000009 GELG00000000012 GELG00000000013
run_mutect GELS00000000010 GELG00000000014 GELG00000000015
