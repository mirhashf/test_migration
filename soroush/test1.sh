#! /bin/bash
readbedfile="true"
while true; do
  read -r lineA <&3
  if [ $readbedfile == "true" ]; then
   read -r lineB <&4
  fi
  if [ -z "$lineA" -o -z "$lineB" ]; then
    break
  fi
  if [[ $lineA == @chr* ]]; then
    readbedfile="true"
    fastq=$(echo ${lineA:1})
    bed=$(echo $lineB | awk '{print $4}')
    echo "1: "$fastq > tmp.txt
    echo "2: "$bed > tmp.txt
    if [[ $fastq == $bed ]]; then
     newcoord=$(echo $lineB | awk '{print ($2+1)"_"($3+1)}')
     concat=$(echo $fastq | awk -v nc=$newcoord '{split($fastq,a,"_");print a[1]"_"a[2]"_"nc}')
     echo $concat >> fastq2.out
    else
     echo $(echo "UNMAPPED_READ"| awk -v $fs '{print "_"fs}') >> fastq2.out
    fi
  else
    echo $lineA >> fastq2.out
    readbedfile="false"
  fi
done 3<fastqoriginal.fastq 4<simulated.lane0.bwa.read1.fastq.LIFTED4M.bed
