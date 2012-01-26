#!/bin/bash
./merge_bam.py s3://seqalto-results/A804NLABXX.s_1_1.fq.gz-A804NLABXX.s_1_2.fq.gz/bwa_sorted.bam \
s3://seqalto-resutls/A804NLABXX.s_2_1.fq.gz-A804NLABXX.s_2_2.fq.gz/bwa_sorted.bam \
s3://seqalto-resutls/A804NLABXX.s_3_1.fq.gz-A804NLABXX.s_3_2.fq.gz/bwa_sorted.bam \
s3://seqalto-resutls/A804NLABXX.s_4_1.fq.gz-A804NLABXX.s_4_2.fq.gz/bwa_sorted.bam \
s3://seqalto-resutls/A804NLABXX.s_5_1.fq.gz-A804NLABXX.s_5_2.fq.gz/bwa_sorted.bam 
