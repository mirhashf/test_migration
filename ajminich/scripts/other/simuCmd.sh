
REF=/mnt/scratch0/public/data/annotation/hg19/chr21.fa

###########
echo "*** simulating reads"
~jianl/work/simuReads/simuRead random_genome $REF /mnt/scratch0/public/data/variants/dbSNP/CEU/CEU-1409-21.vcf /mnt/scratch0/jianl/test/xx 10 0.15
~jianl/work/simuReads/simuRead pair_reads /mnt/scratch0/jianl/test/xx 100 100 7200000 330 7 3 tt

echo "*** bwa align ***"
bwa aln $REF tt_1.fq -q 20 -t 16 1>read1.sai
bwa aln $REF tt_2.fq -q 20 -t 16 1>read2.sai

echo "*** bwa pe ***"
bwa sampe -P $REF read1.sai read2.sai tt_1.fq tt_2.fq > bwaAlign.sam
echo "*** mousetail align ***"

mousetail27 align $REF\_22.midx -m 330 -i 400 --trim 20 -p 2 -1 tt_1.fq -2 tt_2.fq > seqAlign.sam

grep "pa;" seqAlign.sam > paAligned.sam
grep "ma;" seqAlign.sam > maAligned.sam

mousetail27 check_pair /mnt/scratch0/jianl/test/xx.pa.map paAligned.sam 21 3600000 1 1 checkPa 1>seqPa.out 2>seqPa.err
mousetail27 check_pair /mnt/scratch0/jianl/test/xx.ma.map maAligned.sam 21 3600000 1 1 checkMa 1>seqMa.out 2>seqMa.err
#echo "###### seq var#######"

#echo "######bwa var#######"

