#!/bin/bash

set -e

if [ "$(id -u)" != "0" ]; then
    echo "This script would run better as root."
    exit 1
fi

OUTFILEPATH="/var/log/disk-usage"
OUTFILENAME="disk-usage-report_$(date +%Y%m%d)"
OUTFILE=${OUTFILEPATH}/${OUTFILENAME}

touch ${OUTFILE}

SEARCH_DIRS='
/net/kodiak/volumes/lake/shared
/net/kodiak/volumes/lake/shared/users
/net/kodiak/volumes/river/shared
/net/kodiak/volumes/river/shared/users
/net/kodiak/volumes/delta/shared
/net/kodiak/volumes/delta/shared/prj
/net/kodiak/volumes/delta/shared/home
/net/hippo/volumes/wadi/shared
/net/hippo/volumes/wadi/shared/prj'

for DIR in ${SEARCH_DIRS}; do
    pushd ${DIR} 1> /dev/null
    echo "Reporting directory ${DIR}" | tee -a ${OUTFILE}
    du -h --max-depth=1 | sort -hr | tee -a ${OUTFILE}
    popd 1> /dev/null
done

