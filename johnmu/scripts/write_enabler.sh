#!/bin/bash

# Adds group write permissions to a file or directory
# It is careful not to add execute to non-dirctories

if [ -z "$1" ]; then
    echo "Usage: ${0} <file or dir>"
    exit
fi

chmod -R g+wr ${1}
find ${1} -type d -exec chmod g+x {} \;
