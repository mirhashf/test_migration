#!/bin/bash

# Sorts a provided BAM file and extracts the sorted regions.

if [[ ${#} -lt 3 ]]; then
  echo "Region Extraction Script"
  echo "Written by Jian Li and AJ Minich, April-May 2012"
  echo ""
  echo "Use: sudo ${0} <inputBAM> <region> <finalFile>"
  echo "  <inputFile> - raw alignment file (as SAM or BAM) to extract regions from"
  echo "  <region>    - the region to extract"
  echo "  <finalFile> - the final BAM file to write extracted regions to"
  echo ""
  echo "Notes:"
  echo "- Input file will not be sorted if it is already in coordinate order."
  echo "- If sorting is needed, sorted BAM file will be written to <inputFile>.sorted"
  
  exit 0
fi

INPUT_FILE=${1}
REGION=${2}
FINAL_BAM=${3}

PICARD=/mnt/scratch2/public/programs/picard
TMP=/tmp

header=`samtools view -h ${INPUT_FILE}`
if [[ "${header}" == *SO:coordinate* ]]
then
  echo "Input file is already sorted by coordinate.";
  
  SORTED_BAM=${INPUT_FILE}
  
else

  SORTED_BAM=${INPUT_FILE}.sorted

    java -Xms5g -Xmx5g -jar ${PICARD}/dist/SortSam.jar \
        INPUT=${INPUT_FILE} \
        VALIDATION_STRINGENCY=LENIENT \
        OUTPUT=${SORTED_BAM} \
        SORT_ORDER=coordinate \
        TMP_DIR=${TMP}
fi

java -Xms5g -Xmx5g -jar ${PICARD}/dist/BuildBamIndex.jar \
	INPUT=${SORTED_BAM} \
	VALIDATION_STRINGENCY=LENIENT

# Get just the reads aligned to the selected region
echo "Getting reads in the region: ${REGION}"
samtools view -bh ${SORTED_BAM} \
    ${REGION} > ${FINAL_BAM}
