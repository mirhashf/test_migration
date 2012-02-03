# compSNPvcf.R
# VCF File Comparison
# Written by Jian Li
#
# Run as:
#
#   R --slave --args inputFile1 inputFile2 editfile<compSNPvcf.R
#
# - input and editfile should be in VCF format.
# - editfile could be "/mnt/scratch0/public/data/variants/dbSNP/CEU/CEU-1409-21.vcf"
# - outputs three files:
#     > inputFile1 specific prediction
#     > inputFile2 specific SNPs
#     > intersection of the total sets of predictions
#
# Output file: those lines with an "rs" Id in the "ID" column mean they are from dbSNP.

#editfn="/home/public/data/variants/dbSNP/21-1409-CEU.vcf"

Args <- commandArgs()
f1n=Args[4]
f2n=Args[5]

editfn=Args[6]
orig=read.table(editfn,sep="\t",stringsAsFactors=FALSE)
names(orig)=c("CHROM","POS","ID","REF","dbALT","dbQUAL","dbFILTER","dbINFO")

snpFlag=sapply(1:nrow(orig),function(x) grepl("VC=SN",orig[x,8]))# | grepl("VC=MIX",orig[x,8]))
snpInd=which(snpFlag == TRUE)
idInd=which(snpFlag == FALSE)


orig=unique(orig[snpInd,c(1:2)])

f10=read.table(f1n,sep="\t",stringsAsFactors=FALSE)
f20=read.table(f2n,sep="\t",stringsAsFactors=FALSE)


names(f10)=c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","SAMPLE")
names(f20)=c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","SAMPLE")

f1=f10[which(f10$FILTER == "PASS"),]#f10$QUAL >100),]#f1[which(f10$FILTER == "PASS" & gqFlag1 == TRUE & hetFlag1 == TRUE),]
f2=f20[which(f20$FILTER == "PASS"),]#f20$QUAL >100),]#f2[which(f20$FILTER == "PASS" & gqFlag2 == TRUE & hetFlag2 == TRUE),]


f12_both=merge(f1,f2,by=c("CHROM","POS","ID","REF"))
commPos=f12_both$POS
f1_onlyPos=setdiff(f1$POS,commPos)
f2_onlyPos=setdiff(f2$POS,commPos)
f1_only=f1[which(f1$POS %in% f1_onlyPos),] 
f2_only=f2[which(f2$POS	%in% f2_onlyPos),]

write.table(f1_only,file="specSNP1.txt",sep="\t",quote=FALSE,row.names=FALSE)
write.table(f2_only,file="specSNP2.txt",sep="\t",quote=FALSE,row.names=FALSE)


comm=length(commPos)


uniq1T=merge(orig,f1_only,by=c("CHROM","POS"))
uniq2T=merge(orig,f2_only,by=c("CHROM","POS"))

uniq1FP=f1[which(f1$POS %in% (setdiff(f1_onlyPos,uniq1T$POS))),]
uniq2FP=f2[which(f2$POS %in% (setdiff(f2_onlyPos,uniq2T$POS))),]

write.table(f12_both,file="snpInAll.txt",sep="\t",quote=FALSE,row.names=FALSE)
#write.table(uniq1FP,file="seqOnlyFP.txt",sep="\t",quote=FALSE,row.names=FALSE)
#write.table(uniq2FP,file="bwaOnlyFP.txt",sep="\t",quote=FALSE,row.names=FALSE)

fbT=merge(orig,f12_both,by=c("CHROM","POS"))


TP1=dim(fbT)[1]+dim(uniq1T)[1]
TP2=dim(fbT)[1]+dim(uniq2T)[1]
FP1=length(f1_onlyPos)-dim(uniq1T)[1]+comm-dim(fbT)[1]
FP2=length(f2_onlyPos)-dim(uniq2T)[1]+comm-dim(fbT)[1]
FN1=dim(orig)[1]-TP1
FN2=dim(orig)[1]-TP2

PPV1 = TP1 / (TP1 + FP1)
FDR1 = FP1 / (FP1 + TP1)

PPV2 = TP2 / (TP2 + FP2)
FDR2 = FP2 / (FP2 + TP2)

cat("\n")

cat(TP1,paste(TP1/dim(orig)[1]*100,"%"),FP1,FN1,sep=",")
cat("\n")
cat(TP2,paste(TP2/dim(orig)[1]*100,"%"),FP2,FN2,sep=",")
cat("\n")

