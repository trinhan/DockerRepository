#!/usr/bin/Rscript
"usage: \n FindCosmic.R [--maffile=<file> --outputname=<string> --gnomadcutoff=<float> --onlyCoding=<boolean>]
\n options:\n --maffile=<file> maffile that has been annotated using vep & oncokb.
\n --outputname=<string> output string
\n --gnomadcutoff=<float> [default: 0.1]
\n --onlyCoding=<boolean> [default: TRUE]" -> doc

library("docopt", quietly = T)
opts <- docopt(doc)

suppressMessages(library(data.table, quietly = T))
suppressMessages(library(dplyr, quietly = T))
print('RUN COSMIC FILTER')
InputData=fread(opts$maffile, sep="\t")

print('Finding cancer variants')
## Known Variants associated with cancer?
# Find variants with Cancer Gene Census Tier in 1:2
bx1=which(InputData$CancerGeneCensus.Tier%in%c("Hallmark 1", "1", "Hallmark 2", "2")) 
# Find variants with Cancer Gene Census Mut Tier in 1:3
cx1 = which(InputData$CancerMutationCensus.Tier%in%c(1:3))
sprintf('%s in Cosmic Gene Census', length(bx1))
sprintf('%s in Cosmic Mut Census', length(cx1))
dx1=grep("Oncogenic", InputData$Oncokb.ONCOGENIC)
sprintf('%s in Oncokb', length(dx1))
mx1=unique(sort(c(bx1, cx1, dx1)))
sprintf('Total %s cancer variants', length(mx1))
Keep1=InputData[mx1, ]

sprintf('Filter based on gnomad cutoff of %s', opts$gnomadcutoff)
select1=which(Keep1$AF_max<=opts$gnomadcutoff)
Keep1=Keep1[select1, ]
sprintf("%s variants after gnomad filter" ,nrow(Keep1))

if (opts$onlyCoding){
  print('Keep only coding mutations')
  select1=which(Keep1$ConsB==1)
  Keep1=Keep1[select1, ]
  sprintf("%s variants after functional filter" ,nrow(Keep1))
}

Keep1=Keep1%>%select(!c("AF_max"))
print('FINISHED COSMIC SEARCH')
write.table(Keep1, file=paste(opts$outputname, "cancerVariants.filt.maf", sep=""), sep="\t", row.names = F, quote = F)

