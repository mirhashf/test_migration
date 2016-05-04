#!/bin/bash

function hashname {
    echo ${1} | sha1sum | cut -c 1-8
}

for S in ${@}; do
    pushd ${S}/Fastq

    for F in $(find . -maxdepth 1 -type f -name "*.fastq.gz"); do
        ~/s3cmd/s3cmd put ${F} s3://bina.va.uswest1.inputdata/$(hashname ${S})/$(echo ${F} | cut -d"/" -f2-)
    done

    popd

    shift
done

