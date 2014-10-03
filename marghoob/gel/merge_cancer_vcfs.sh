#!/bin/bash

#500490e8-24c1-402a-877a-b79fa5c8a14b
#0e52a8fa-a3d8-4239-aa66-fc47d6832625
#f69a1798-8d16-48ca-90ea-d15056f5d30a
#6b5071c6-6375-4e3d-9641-d63b986d41bc

TABIX=~/lake/opt/tabix/tabix
BGZIP=~/lake/opt/tabix/bgzip
merge_vcf=../merge_vcfs.sh
for jobid in 500490e8-24c1-402a-877a-b79fa5c8a14b 0e52a8fa-a3d8-4239-aa66-fc47d6832625 f69a1798-8d16-48ca-90ea-d15056f5d30a 6b5071c6-6375-4e3d-9641-d63b986d41bc; do
  jobdir=/net/kodiak/volumes/delta/shared/prj/GELA20140221/cancer/runs/mutect/GELS000000000??/$jobid
  jobdir=$(echo $jobdir)
  for vartype in snp indel; do
    ls -l $jobdir/vcf/$vartype/ALL.vcf.gz
    vcf-validator $jobdir/vcf/$vartype/ALL.vcf.gz &>$jobid.$vartype.validator.log &
    continue
    VCFS_DIR=$jobdir/vcf/$vartype
    REFERENCE=~/lake/references/human/hg19/hg19.fa
    [ -w "$VCFS_DIR" ] && echo "Can write to $VCFS_DIR"
    [ -e "$VCFS_DIR/1.vcf.gz" ] && REFERENCE=~/lake/references/human/human_g1k_v37_decoy/human_g1k_v37_decoy.fasta
    VCFS_DIR=$VCFS_DIR TABIX=$TABIX BGZIP=$BGZIP REFERENCE=$REFERENCE OUTPUT_VCF=$jobdir/vcf/$vartype/ALL.vcf $merge_vcf
  done
done
wait
