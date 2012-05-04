#!/bin/bash

# Runs Snyder on the remaining lanes.

for i in {1..5}
do

  ~/sandbox/jianl/seqaltoAlign.ec2.sh A808HKABXX.s_${i} SLB

  s3cmd put ~/execution/seqalto_A808HKABXX.s_${i}.chr1_AG.bam s3://bina.results/snyder/
  s3cmd put ~/execution/seqalto_A808HKABXX.s_${i}.log s3://bina.results.snyder/

  # Cleanup
  rm -f ~/execution/*.sam
  rm -f ~/execution/*${i}.bam
  rm -f ~/data/*.fq

done

