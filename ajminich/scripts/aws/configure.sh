#! /bin/bash

# Setup script for running SeqAlto on EC2
# Note: run this script with sudo so it can create and write directories.

# Get the most recent SeqAlto binary and set up the symbolic link
echo "Setting up most recent version of SeqAlto aligner from S3."
s3cmd get s3://bina.programs/aligner.aws ~/programs/aligner.aws
chmod a+x ~/programs/aligner.aws
rm ~/bin/seqalto
ln -s ~/bin/seqalto ~/programs/aligner.aws

# Set up /ebs
echo "Setting up EBS persistent storage."
chmod a+rw /ebs
mkdir /ebs/data
mkdir /ebs/execution

# Set up /mnt
chmod a+rw /mnt

# Set up needed files
echo "Getting FASTA Index and generating Seqalto index."
s3cmd get s3://bina.data/hg19/ucsc.hg19.fasta /ebs/data/hg19.fa
seqalto -mode index \
    -ref /ebs/data/hg19.fa \
    -index_mode 0 \
    -kmer_size 22 \
    -index_name /ebs/data/hg19.fa_22.sidx \
    -logtostderr






