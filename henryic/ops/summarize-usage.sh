#!/bin/bash

if [ "$(id -u)" != "0" ]; then
    echo "This script would run better as root."
    exit 1
fi

# Show overall river usage breakdown
pushd /net/kodiak/volumes/river/shared 1> /dev/null
echo "river:"
du -h --max-depth=1
popd 1> /dev/null

# Show 15 largest directories in river/output
pushd /net/kodiak/volumes/river/shared/output 1> /dev/null
echo "river/output:"
echo "Top 15 job output directories:"
OUTPUT_USAGE=$(mktemp)
du -h --max-depth 1 | sort -rh | head -n 16 | tail -n 15 > ${OUTPUT_USAGE}
cat ${OUTPUT_USAGE}
cat ${OUTPUT_USAGE} | awk '{ print $2 }' | xargs ls -ld
rm ${OUTPUT_USAGE}
unset OUTPUT_USAGE
popd 1> /dev/null

# Show overall river user directory usage
pushd /net/kodiak/volumes/river/shared/users 1> /dev/null
echo "river/users:"
du -h --max-depth 1 | sort -rh
popd 1> /dev/null

# Show overall lake usage breakdown
pushd /net/kodiak/volumes/lake/shared 1> /dev/null
echo "lake:"
du -h --max-depth=1
popd 1> /dev/null

# Show overall lake user directory usage
pushd /net/kodiak/volumes/lake/shared/users 1> /dev/null
echo "lake/users:"
du -h --max-depth 1 | sort -rh
popd 1> /dev/null

