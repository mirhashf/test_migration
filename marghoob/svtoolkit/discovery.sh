#!/bin/bash

# If you adapt this script for your own use, you will need to set these two variables based on your environment.
# SV_DIR is the installation directory for SVToolkit - it must be an exported environment variable.
# SV_TMPDIR is a directory for writing temp files, which may be large if you have a large data set.
export LAKE=/net/kodiak/volumes/lake/shared
export SV_DIR=$LAKE/opt/svtoolkit/
SV_TMPDIR=./tmpdir
reference=$LAKE/references/human/human_g1k_v37_decoy/human_g1k_v37_decoy.fasta
genome_mask=/net/kodiak/volumes/delta/shared/ftp/ftp.broadinstitute.org/pub/svtoolkit/svmasks/Homo_sapiens_assembly19.mask.101.fasta
ploidy_map=/net/kodiak/volumes/delta/shared/ftp/ftp.broadinstitute.org/pub/svtoolkit/ploidymaps/humgen_g1k_v37_ploidy.map
copy_number_mask=/net/kodiak/volumes/delta/shared/ftp/ftp.broadinstitute.org/pub/svtoolkit/cn2masks/cn2_mask_g1k_v37.fasta
gender_map=/net/kodiak/volumes/river/shared/users/marghoob/sabeti/sry/gender.map
genstrip_conf=$SV_DIR/conf/genstrip_parameters.txt
regions="-L 22"

runDir=$1
shift

bams_opt=
for i in $*; do
  bams_opt="$bams_opt -I $i"
done

[ -z "$runDir" -o -z "$bams_opt" ] && echo "./discovery.sh <runDir> <bams>" && exit 1

mkdir -pv $runDir

sites=$runDir/discovery.vcf
genotypes=$runDir/genotypes.vcf

# These executables must be on your path.
which java > /dev/null || exit 1
which Rscript > /dev/null || exit 1
which samtools > /dev/null || exit 1

# For SVAltAlign, you must use the version of bwa compatible with Genome STRiP.
export PATH=${SV_DIR}/bwa:${PATH}
export LD_LIBRARY_PATH=${SV_DIR}/bwa:${LD_LIBRARY_PATH}

mx="-Xmx4g"
classpath="${SV_DIR}/lib/SVToolkit.jar:${SV_DIR}/lib/gatk/GenomeAnalysisTK.jar:${SV_DIR}/lib/gatk/Queue.jar"

mkdir -p ${runDir}/logs || exit 1
mkdir -p ${runDir}/metadata || exit 1

# Display version information.
java -cp ${classpath} ${mx} -jar ${SV_DIR}/lib/SVToolkit.jar

# Run preprocessing.
# For large scale use, you should use -reduceInsertSizeDistributions, but this is too slow for the installation test.
# The method employed by -computeGCProfiles requires a CN2 copy number mask and is currently only supported for human genomes.
java -cp ${classpath} ${mx} \
    org.broadinstitute.sting.queue.QCommandLine \
    -S ${SV_DIR}/qscript/SVPreprocess.q \
    -S ${SV_DIR}/qscript/SVQScript.q \
    -gatk ${SV_DIR}/lib/gatk/GenomeAnalysisTK.jar \
    --disableJobReport \
    -cp ${classpath} \
    -configFile $genstrip_conf \
    -tempDir ${SV_TMPDIR} \
    -R $reference \
    -genomeMaskFile $genome_mask \
    -ploidyMapFile $ploidy_map \
    -copyNumberMaskFile $copy_number_mask \
    -genderMapFile $gender_map \
    -runDirectory ${runDir} \
    -md ${runDir}/metadata \
    -disableGATKTraversal \
    -useMultiStep \
    -computeGCProfiles \
    -jobLogDir ${runDir}/logs \
    $bams_opt \
    $regions \
    -run \
    -useMultiStep \
    || exit 1

if ((0)); then
java -cp ${classpath} ${mx} \
    org.broadinstitute.sting.queue.QCommandLine \
    -S ${SV_DIR}/qscript/SVAltAlign.q \
    -S ${SV_DIR}/qscript/SVQScript.q \
    -gatk ${SV_DIR}/lib/gatk/GenomeAnalysisTK.jar \
    -configFile $genstrip_conf \
    -tempDir ${SV_TMPDIR} \
    -md ${runDir}/metadata \
    -R $reference \
    -vcf input.vcf \
    $bams_opt \
    -O ${runDir}/alt.bam \
    -run \
    $regions \
    -jobLogDir ${runDir}/logs || exit 1
fi

# Run discovery.
java -cp ${classpath} ${mx} \
    org.broadinstitute.sting.queue.QCommandLine \
    -S ${SV_DIR}/qscript/SVDiscovery.q \
    -S ${SV_DIR}/qscript/SVQScript.q \
    -gatk ${SV_DIR}/lib/gatk/GenomeAnalysisTK.jar \
    --disableJobReport \
    -cp ${classpath} \
    -configFile $genstrip_conf \
    -tempDir ${SV_TMPDIR} \
    -R $reference \
    -genomeMaskFile $genome_mask \
    -genderMapFile $gender_map \
    -runDirectory ${runDir} \
    -md ${runDir}/metadata \
    -disableGATKTraversal \
    -jobLogDir ${runDir}/logs \
    $regions \
    -minimumSize 100 \
    -maximumSize 1000000 \
    -suppressVCFCommandLines \
    $bams_opt \
    -O ${sites} \
    -run \
    || exit 1

# Run genotyping on the discovered sites.
java -cp ${classpath} ${mx} \
    org.broadinstitute.sting.queue.QCommandLine \
    -S ${SV_DIR}/qscript/SVGenotyper.q \
    -S ${SV_DIR}/qscript/SVQScript.q \
    -gatk ${SV_DIR}/lib/gatk/GenomeAnalysisTK.jar \
    --disableJobReport \
    -cp ${classpath} \
    -configFile $genstrip_conf \
    -tempDir ${SV_TMPDIR} \
    -R $reference \
    -genomeMaskFile $genome_mask \
    -genderMapFile $gender_map \
    -runDirectory ${runDir} \
    -md ${runDir}/metadata \
    -disableGATKTraversal \
    -jobLogDir ${runDir}/logs \
    $bams_opt \
    -vcf ${sites} \
    -O ${genotypes} \
    $regions \
    -run \
    || exit 1

