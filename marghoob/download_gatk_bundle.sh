#!/bin/bash

FTP_BASE="ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/2.5/"

bundleid=$1
outdir=$2

#wget -P $outdir $FTP_BASE/$bundleid/*.fasta.*
#wget -P $outdir $FTP_BASE/$bundleid/*.dict.*

for prefix in 1000G dbsnp hapmap Mills
do
wget -P $outdir $FTP_BASE/$bundleid/$prefix*
done

gunzip $outdir/*.gz 
