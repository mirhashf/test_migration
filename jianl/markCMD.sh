
TMP=./
PICARD=~/downloads/picard/dist
d=$1
f=$2

java -Xms5g -Xmx5g -jar $PICARD/BuildBamIndex.jar INPUT=$d/$f.bam VALIDATION_STRINGENCY=LENIENT

echo ">>> Marking duplicates"

java -Xms5g -Xmx5g -jar $PICARD/MarkDuplicates.jar \
        TMP_DIR=$TMP \
        I=$d/$f.bam\
        O=$d/$f\_marked.bam\
        M=$d/$f.metrics \
        VALIDATION_STRINGENCY=SILENT \
        ASSUME_SORTED=true \
        REMOVE_DUPLICATES=true \
        READ_NAME_REGEX=".+:[0-9]+:([0-9]+):([0-9]+):([0-9]+).*" 

#        OPTICAL_DUPLICATE_PIXEL_DISTANCE=100

echo "*** Finished removing duplicates ***"

#java -Xms5g -Xmx5g -jar $PICARD/BuildBamIndex.jar INPUT=$d/$f\_marked.bam VALIDATION_STRINGENCY=LENIENT

echo "DONE"

