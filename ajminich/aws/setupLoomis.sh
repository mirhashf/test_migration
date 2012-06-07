# !/bin/bash

# Setup script for running Loomis on EC2.
# Should be run as root.

# For more information, please visit the Loomis guide at
# https://binatechnologies.atlassian.net/wiki/display/SEQ/Loomis+User%27s+Guide

# Directories and files
HOME=/home/ubuntu/
ZOOKEEPER=${HOME}/zookeeper-3.3.5/
SEQALTO=${HOME}/seqalto/
DATA=/mnt/data/
EXECUTION=/mnt/execution/

# Indexing options
INDEX_MODE=0
KMER_SIZE=22

# Build necessary files
cd ${SEQALTO}
scons aligner -j24
scons samsorter -j24
ant dist -file loomis/server/build.xml

# Get data
mkdir ${DATA}
chmod a+rw ${DATA}
cd ${DATA}
s3cmd get s3://bina.data/hg19/ucsc.hg19.fasta
samtools faidx ucsc.hg19.fasta

s3cmd get s3://bina.data/dbSNP/dbsnp.hg19.vcf
s3cmd get s3://bina.data/svdata/windows.tar.gz
tar -xf windows.tar.gz

${SEQALTO}/aligner/build/aligner -mode index \
    -ref ${DATA}/ucsc.hg19.fasta
    -index_mode ${INDEX_MODE} \
    -kmer_size ${KMER_SIZE} \
    -index_name ${DATA}/ucsc.hg19.fasta_${KMER_SIZE}.midx \
    -logtostderr=1

# Start Zookeeper
${ZOOKEEPER}/bin/zkServer.sh start

# Setup complete
echo "Setup complete: all resources successfully downloaded/initialized."