#!/bin/bash

fileList=`ls $1*.sam`
outFile=$2

catArgs=""
for file in $fileList
do
	/home/ubuntu/programs/convert_snap_sam.pl $file > $file.fixed.sam
	catArgs="${catArgs} I=${file}.fixed.sam"
done

java -Djava.io.tmpdir="/mnt/execution/.tmp/" -Xms2g -Xmx4g -jar /home/ubuntu/programs/picard/dist/MergeSamFiles.jar $catArgs O=$outFile

rm $1*.sam.fixed.sam
