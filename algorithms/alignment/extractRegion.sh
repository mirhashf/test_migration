#!/bin/bash

# Sorts a provided BAM file and extracts the sorted regions.

if [[ ${#} -lt 4 ]]; then
  echo "Region Extraction Script"
  echo "Written by Jian Li and AJ Minich, April-May 2012"
  echo ""
  echo "Use: sudo ${0} <inputBAM> <sortedBAM> <region> <finalFile>"
  echo "  <inputFile> - raw alignment file (as SAM or BAM) to "
  
  exit 0
fi

INPUT_FILE=${1}
SORTED_BAM=${2}
REGION=${3}
FINAL_BAM=${4}

header=`/bin/samtools view -h ${INPUT_FILE}`
if [[ "${header}" == *SO:coordinate* ]]
then
  echo "Input file is already sorted by coordinate.";
else

    sudo java -Xms5g -Xmx5g -jar ~/programs/picard/dist/SortSam.jar \
        INPUT=${INPUT_FILE} \
        VALIDATION_STRINGENCY=LENIENT \
        OUTPUT=${SORTED_BAM} \
        SORT_ORDER=coordinate \
        TMP_DIR=$TMP
fi

sudo java -Xms5g -Xmx5g -jar ~/programs/picard/dist/BuildBamIndex.jar \
	INPUT=${SORTED_BAM} \
	VALIDATION_STRINGENCY=LENIENT

# Get just the reads aligned to the selected region
echo "Getting reads in the region: ${REGION}"
sudo /bin/samtools view -bh ${SORTED_BAM} \
    ${REGION} > ${FINAL_BAM}
