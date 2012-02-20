#! /bin/bash

# configure.sh
# Configures the scratchdisk.

DATA_FOLDER="/mnt/data"
mkdir -p ${DATA_FOLDER}
chmod a+rw ${DATA_FOLDER}
ln -s ${DATA_FOLDER} /home/ubuntu/

EXEC_FOLDER="/mnt/execution"
mkdir ${EXEC_FOLDER}
chmod a+rw ${EXEC_FOLDER}
ln -s ${EXEC_FOLDER}  /home/ubuntu/
