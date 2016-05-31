#!/bin/bash

[ ! -z ${1} ] && KEEP_DAYS=${1} || KEEP_DAYS=90

pushd /net/kodiak/volumes/river/shared/.bina

for f in $(sudo find . -maxdepth 2 -mtime +${KEEP_DAYS} -iname job.json); do if sudo grep -l "\"uri\" *: *\"mount:river\/qe\"" ${f}; then sudo mv $(echo ${f} | cut -d'/' -f2) ../.bina.qe.old/. ; fi; done

popd

