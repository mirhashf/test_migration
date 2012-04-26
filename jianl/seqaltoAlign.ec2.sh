# SeqAlto Alignment for Snyder Reads Data

if [[ ${#} -lt 2 ]]; then
  echo "SeqAlto Alignment Script"
  echo "Written by Jian Li, April 2012"
  echo ""
  echo "Use: sudo ${0} <read_prefix>_#.fq.gz <library>"
  echo "  <read_prefix> - the read from s3://seqalto/ to use (for example, A804NLABXX.s_5)"
  echo "  <library>     - the library to use in add/replace groups (SLB for saliva, BLB for blood)"

  return
fi

TMP=~/mnt

# Get the reads and trim them to quality level 30
sudo s3cmd get s3://seqalto/$1_*.fq.gz ~/data/
sudo gunzip ~/data/$1_*.fq.gz
sudo ~/bin/trim ~/data/$1 30 

# Run the alignment
sudo time ~/bin/seqalto -mode align \
	-1 ~/data/$1_trimmed_1.fq \
	-2 ~/data/$1_trimmed_2.fq \
	-idx ~/data/hg19.fa_22.sidx \
	--template_len_comp_method 2 \
	--enable_batch_pairing \
	-p 8 \
	-h -1 -nw_disable_match_at_ends \
	-verbose \
	2>&1 \
	1>/ebs/execution/seqalto_$1.sam | tee seqalto_$1.log

# Convert to BAM and sort by coordinate
sudo /bin/samtools view -bS /ebs/execution/seqalto_$1.sam > /ebs/execution/seqalto_$1.bam
sudo java -Xms5g -Xmx5g -jar ~/programs/picard/dist/SortSam.jar \
	INPUT=/ebs/execution/seqalto_$1.bam \
	VALIDATION_STRINGENCY=LENIENT \
	OUTPUT=/ebs/execution/seqalto_$1.csorted.bam \
	SORT_ORDER=coordinate

sudo java -Xms5g -Xmx5g -jar ~/programs/picard/dist/BuildBamIndex.jar \
	INPUT=/ebs/execution/seqalto_$1.csorted.bam \
	VALIDATION_STRINGENCY=LENIENT

# Get just the reads aligned to chr1
sudo /bin/samtools view -h /ebs/execution/seqalto_$1.csorted.bam chr1 > /ebs/execution/seqalto_$1.chr1.sam
sudo /bin/samtools view -bS /ebs/execution/seqalto_$1.chr1.sam > /ebs/execution/seqalto_$1.chr1.bam

# Perform add/replace read groups
java -Xms5g -Xmx5g -jar ~/programs/picard/dist/AddOrReplaceReadGroups.jar \
	I=/ebs/execution/seqalto_$1.chr1.bam \
	O=/ebs/execution/seqalto_$1.chr1_AG.bam  \
	SORT_ORDER=coordinate \
	RGID=$1 \
	RGLB=$2 \
	RGPL=ILLUMINA \
	RGSM=$2 \
	RGPU=HiSeq \
	VALIDATION_STRINGENCY=LENIENT  \
	TMP_DIR=$TMP

echo "Execution complete."

