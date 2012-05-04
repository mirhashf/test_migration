#!/bin/bash
# SeqAlto Alignment for Snyder Reads Data

if [[ ${#} -lt 2 ]]; then
  echo "SeqAlto Alignment Script"
  echo "Written by Jian Li, April 2012"
  echo "Modified by AJ Minich, May 2012"
  echo ""
  echo "Use: sudo ${0} <read_prefix>_#.fq.gz <library>"
  echo "  <read_prefix> - the read from s3://seqalto/ to use (for example, A804NLABXX.s_5)"
  echo "  <library>     - the library to use in add/replace groups (SLB for saliva, BLB for blood)"

  return
fi

LANE=${1}
LIB=${2}

TMP=/mnt
QUAL=30
REGION=chr1
DATA_FOLDER=/snyder/data/
EXECUTION_FOLDER=.

SEQALTO_OUT_SAM=${EXECUTION_FOLDER}/seqalto_${LANE}.sam
SEQALTO_OUT_BAM=${EXECUTION_FOLDER}/seqalto_${LANE}.bam
SORTED_BAM=${EXECUTION_FOLDER}/seqalto_${LANE}.sorted.bam
SORTED_REGION_BAM=${EXECUTION_FOLDER}/seqalto_${LANE}.${REGION}.csorted.bam
SORTED_AG_REGION_BAM=${EXECUTION_FOLDER}/seqalto_${LANE}.${REGION}_AG.bam

echo "Running Seqalto alignment on Snyder lane ${LANE} from sample {$LIB}."

# Get the reads and trim them to quality level 30
s3cmd get s3://seqalto/${LANE}_*.fq.gz ${DATA_FOLDER}

echo "Unzipping lanes..."
gunzip ${DATA_FOLDER}/${LANE}_*.fq.gz

echo "Trimming lanes..."
~/sandbox/algorithms/alignment/trim.sh ${DATA_FOLDER}/$1 ${QUAL}

# Remove the original FASTQ files
rm ${DATA_FOLDER}/${LANE}_1.fq
rm ${DATA_FOLDER}/${LANE}_2.fq

# Run the alignment
sudo time ~/bin/seqalto -mode align \
	-1 ${DATA_FOLDER}/${LANE}_trimmed_1.fq \
	-2 ${DATA_FOLDER}/${LANE}_trimmed_2.fq \
	-idx /ebs/data/hg19.major.fa_22.sidx \
	--template_len_comp_method 2 \
	--enable_batch_pairing \
	-p 8 \
	-h -1 -nw_disable_match_at_ends \
	-verbose \
	2>&1 \
	1>${SEQALTO_OUT_SAM} | tee ${EXECUTION_FOLDER}/seqalto_${LANE}.log

# Convert to BAM and sort by coordinate
echo "Converting to BAM..."
sudo /bin/samtools view -bS ${SEQALTO_OUT_SAM} > ${SEQALTO_OUT_BAM}
rm ${SEQALTO_OUT}

sudo java -Xms5g -Xmx5g -jar ~/programs/picard/dist/SortSam.jar \
	INPUT=${SEQALTO_OUT_BAM} \
	VALIDATION_STRINGENCY=LENIENT \
	OUTPUT=${SORTED_BAM} \
	SORT_ORDER=coordinate \
	TMP_DIR=$TMP

# Remove now-unnecessary raw BAM file
rm ${SEQALTO_OUT_BAM}

sudo java -Xms5g -Xmx5g -jar ~/programs/picard/dist/BuildBamIndex.jar \
	INPUT=${SORTED_BAM} \
	VALIDATION_STRINGENCY=LENIENT

# Get just the reads aligned to the selected region
echo "Getting reads in the region: ${REGION}"
sudo /bin/samtools view -bh ${SORTED_BAM} \
    ${REGION} > ${SORTED_REGION_BAM}

# Perform add/replace read groups
java -Xms5g -Xmx5g -jar ~/programs/picard/dist/AddOrReplaceReadGroups.jar \
	I=${SORTED_REGION_BAM} \
	O=${SORTED_AG_REGION_BAM} \
	SORT_ORDER=coordinate \
	RGID=${LANE} \
	RGLB=${LIB} \
	RGPL=ILLUMINA \
	RGSM=${LIB} \
	RGPU=HiSeq \
	VALIDATION_STRINGENCY=LENIENT  \
	TMP_DIR=${TMP}

echo "Execution complete."

