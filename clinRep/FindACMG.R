#!/usr/bin/Rscript
"usage: \n FindACMG.R [--maffile=<file> --outputname=<string> --ACMG=<file>]
\n options:\n --maffile=<file> maffile that has been annotated using vep & oncokb.
\n --outputname=<string> output string
\n --ACMG=<file> Table of ACMG variants" -> doc

library("docopt", quietly = T)
opts <- docopt(doc)

suppressMessages(library(data.table, quietly = T))
suppressMessages(library(dplyr, quietly = T))
print('Run ACMG filter')
InputData=fread(opts$maffile, sep="\t")

print('Read in ACMG file')
## ACMG variants
ACMGtable=read.csv(opts$ACMG)
Lx1=which(InputData$SYMBOL%in%ACMGtable$Gene)
strSplitA=strsplit(InputData$ClinVar.Sig, "&")
Nx2=sapply(strSplitA, function(x) length(which(c("pathogenic", "Pathogenic", "Likely_pathogenic", "likely_pathogenic")%in%x)))
Lx2=which(Nx2>0)
Lx3=intersect(Lx1, Lx2)
sprintf("%s variants in ACMG list", length(Lx3))
Keep0=InputData[Lx3, ]%>%select(!"AF_max")

print('Finished! Writing out files')
write.table(Keep0, file=paste(opts$outputname, "ACMG.filt.maf", sep=""), sep="\t", row.names = F, quote = F)
