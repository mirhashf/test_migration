#!/bin/bash

#generate the two binaries for both the codes
run_fail=1
ungapped=0
single_gap=1
paired=1
genome=genome_all

echo "generating the binaries"
make clean
make DEBUG=1 SINGLE_GAP=$single_gap 
mv seqalto seqalto-debug
if [ $run_fail -gt 0 ]
then
  make clean
  make DEBUG=0 SINGLE_GAP=$single_gap
fi

cd ../aligner-basic/changes
make clean 
make DEBUG=1 SINGLE_GAP=$single_gap
mv mousetail mousetail-debug

if [ $run_fail -gt 0 ]
then
  make clean
  make DEBUG=0 SINGLE_GAP=$single_gap
fi

cd ../../aligner-fast

for reads_prefix in test/lane_1_sampled test/real test/reads
do
  align_prefix=$reads_prefix"_align"
  fails=0
  if [ $run_fail -gt 0 ]
  then
    GENOME=$genome READS_PREFIX=$reads_prefix ALIGN_PREFIX=$align_prefix LOG_PREFIX=$reads_prefix DEBUG=0 UNGAPPED=$ungapped SINGLE_GAP=$single_gap PAIRED=$paired make runtest
    diff $align_prefix".sam" $align_prefix"_ref.sam"|grep -e "^<"|grep -v seqalto|cut -f2 -d'<'|cut -f1|sed 's/ //'|uniq > $reads_prefix"_fail.txt"
    echo `wc -l $reads_prefix"_fail.txt"` reads failed
    rm -f $reads_prefix"_fail"*.fq
    for readid in `cat $reads_prefix"_fail.txt"`
    do
      grep -A3 "$readid" $reads_prefix"_1.fq" >> $reads_prefix"_fail_1.fq"
      grep -A3 "$readid" $reads_prefix"_2.fq" >> $reads_prefix"_fail_2.fq"
      let fails=fails+1
      if [ $fails -eq 10 ]
      then 
        break
      fi
    done
    echo $fails failing reads to be tested
  fi 
  if [ $run_fail -eq 0 ]
  then
    fails=1
  fi
  if [ $fails -gt 0 ]
  then
    echo "will run on failing reads now"
    GENOME=$genome READS_PREFIX=$reads_prefix"_fail" ALIGN_PREFIX=$reads_prefix"_fail_align" LOG_PREFIX=$reads_prefix"_fail" DEBUG=1 UNGAPPED=$ungapped SINGLE_GAP=$single_gap PAIRED=$paired make runtest
    echo `diff $reads_prefix"_fail_align.sam" $reads_prefix"_fail_align_ref.sam"|grep -e "^<"|grep -v seqalto|wc -l` differ
  fi
done

