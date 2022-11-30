#!/usr/bin/Rscript
"usage: \n FindDrugAssoc.R [--maffile=<file> --outputname=<string>]
\n options:\n --maffile=<file> maffile that has been annotated using vep & oncokb.
\n --outputname=<string> output string" -> doc

library("docopt", quietly = T)
opts <- docopt(doc)

suppressMessages(library(data.table, quietly = T))
suppressMessages(library(dplyr, quietly = T))
print('RUN CLINVAR DRUG FILTER')
InputData=fread(opts$maffile, sep="\t")

print('Finding drug assoc variants')
## Known Variants associated with drug response or 
nx1=grep("drug_response|protective", InputData$ClinVar.Sig)
sprintf("%s variants related to drug response or protective effect", length(nx1))
Keep2=InputData[nx1, ]%>%select(!"AF_max")
print('FINISHED CLINVAR DRUG FILTER')
write.table(Keep2, file=paste(opts$outputname, "drugprotective.filt.maf", sep=""), sep="\t", row.names = F, quote = F)
