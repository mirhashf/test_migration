#!/bin/bash

for vcfFile in $( find . | grep LP.*\.vcf$ ); do
    tarfile=`echo $vcfFile | sed 's/.*\(LP.*\)\/\(.*\)\/.*/\1\.\2\.tar\.gz/g'`
    tar -czf $tarfile $vcfFile
    s3cmd put $tarfile s3://binabox-redwood/logs/$tarfile
    rm $tarfile
done

