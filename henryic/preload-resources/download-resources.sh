#!/bin/bash

RESOURCE_DIR=resources
DOWNLOAD_DIR=raw
BROAD_BUNDLE_VER=2.8


function download_broad {
    mkdir -pv ftp.broadinstitute.org/bundle/${BROAD_BUNDLE_VER}/${1}
    pushd ftp.broadinstitute.org/bundle/${BROAD_BUNDLE_VER}/${1}
    wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/${BROAD_BUNDLE_VER}/${1}/
    grep "a href" index.html | sed -e 's/.*<a href="\(.*\)">.*/\1/g' | grep -v \.bam | grep -v \.done | wget -N -i -

    for f in $(ls *.md5); do
        echo -n "Verifying file integrity for file $(basename ${f} .md5) ... "
        if [ "$(awk '{ print $1 }' ${f})" != "$(md5sum $(basename ${f} .md5) | awk '{ print $1 }')" ]; then
            echo "MD5SUM mismatch on file ${f}"
        else
            echo "Passed."
        fi
    done

    popd
}

function unpack_broad {
    pushd ftp.broadinstitute.org/bundle/${BROAD_BUNDLE_VER}/${1}
    for f in $(find . -iname "*.vcf.gz" -o -iname "*.fasta.gz" -o -iname "*.fa.gz"); do
        gzip -vdc ${f} > $(basename ${f} .gz)
    done
}

function link_broad {
    pushd ${1}

    find ../raw/ftp.broadinstitute.org/bundle/${BROAD_BUNDLE_VER}/${1}/ -iname "*.fasta" -o -iname "*.fa" -exec ln -s '{}' \;
    find ../raw/ftp.broadinstitute.org/bundle/${BROAD_BUNDLE_VER}/${1}/ -iname "*dbsnp*.vcf" -exec ln -s '{}' \;

    cd vqsr
    find ../../raw/ftp.broadinstitute.org/bundle/${BROAD_BUNDLE_VER}/${1}/ -iname "*.vcf" -exec ln -s '{}' \;
    cd ..

    popd
}


mkdir -pv ${RESOURCE_DIR}/${DOWNLOAD_DIR}
pushd ${RESOURCE_DIR}/${DOWNLOAD_DIR}

# Chromosome list: 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y,MT,GL000207.1,GL000226.1,GL000229.1,GL000231.1,GL000210.1,GL000239.1,GL000235.1,GL000201.1,GL000247.1,GL000245.1,GL000197.1,GL000203.1,GL000246.1,GL000249.1,GL000196.1,GL000248.1,GL000244.1,GL000238.1,GL000202.1,GL000234.1,GL000232.1,GL000206.1,GL000240.1,GL000236.1,GL000241.1,GL000243.1,GL000242.1,GL000230.1,GL000237.1,GL000233.1,GL000204.1,GL000198.1,GL000208.1,GL000191.1,GL000227.1,GL000228.1,GL000214.1,GL000221.1,GL000209.1,GL000218.1,GL000220.1,GL000213.1,GL000211.1,GL000199.1,GL000217.1,GL000216.1,GL000215.1,GL000205.1,GL000219.1,GL000224.1,GL000223.1,GL000195.1,GL000212.1,GL000222.1,GL000200.1,GL000193.1,GL000194.1,GL000225.1,GL000192.1,NC_007605,(hs37d5)
[ "$1" == "--no-download" ] || download_broad b37

# Chromosome list: chrM,chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY,chr1_gl000191_random,chr1_gl000192_random,chr4_ctg9_hap1,chr4_gl000193_random,chr4_gl000194_random,chr6_apd_hap1,chr6_cox_hap2,chr6_dbb_hap3,chr6_mann_hap4,chr6_mcf_hap5,chr6_qbl_hap6,chr6_ssto_hap7,chr7_gl000195_random,chr8_gl000196_random,chr8_gl000197_random,chr9_gl000198_random,chr9_gl000199_random,chr9_gl000200_random,chr9_gl000201_random,chr11_gl000202_random,chr17_ctg5_hap1,chr17_gl000203_random,chr17_gl000204_random,chr17_gl000205_random,chr17_gl000206_random,chr18_gl000207_random,chr19_gl000208_random,chr19_gl000209_random,chr21_gl000210_random,chrUn_gl000211,chrUn_gl000212,chrUn_gl000213,chrUn_gl000214,chrUn_gl000215,chrUn_gl000216,chrUn_gl000217,chrUn_gl000218,chrUn_gl000219,chrUn_gl000220,chrUn_gl000221,chrUn_gl000222,chrUn_gl000223,chrUn_gl000224,chrUn_gl000225,chrUn_gl000226,chrUn_gl000227,chrUn_gl000228,chrUn_gl000229,chrUn_gl000230,chrUn_gl000231,chrUn_gl000232,chrUn_gl000233,chrUn_gl000234,chrUn_gl000235,chrUn_gl000236,chrUn_gl000237,chrUn_gl000238,chrUn_gl000239,chrUn_gl000240,chrUn_gl000241,chrUn_gl000242,chrUn_gl000243,chrUn_gl000244,chrUn_gl000245,chrUn_gl000246,chrUn_gl000247,chrUn_gl000248,chrUn_gl000249
[ "$1" == "--no-download" ] || download_broad hg19

# Chromosome list: 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y,MT
[ "$1" == "--no-download" ] || wget -r ftp://igenome:G3nom3s4u@ussd-ftp.illumina.com/Homo_sapiens/Ensembl/GRCh37/Homo_sapiens_Ensembl_GRCh37.tar.gz

# Chromosome list: chrM,chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY
[ "$1" == "--no-download" ] || wget -r ftp://igenome:G3nom3s4u@ussd-ftp.illumina.com/Homo_sapiens/UCSC/hg19/Homo_sapiens_UCSC_hg19.tar.gz

popd

pushd ${RESOURCE_DIR}

mkdir -pv b37-decoy/bplib       b37-decoy/regions       b37-decoy/somatic       b37-decoy/vqsr
mkdir -pv b37-no-decoy/bplib    b37-no-decoy/regions    b37-no-decoy/somatic    b37-no-decoy/vqsr
mkdir -pv hg19/bplib            hg19/regions            hg19/somatic            hg19/vqsr

mkdir -pv Ensembl_GRCh37
mkdir -pv UCSC_hg19

link_broad hg19
link_broad b37

popd

