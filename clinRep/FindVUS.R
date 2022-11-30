#!/usr/bin/Rscript
"usage: \n FindVUS.R [--maffile=<file> --outputname=<string> --ACMG=<file> --pathwayList=<string> --pathwayFile=<file> --caddscore=<float> --gnomadcutoff=<float> --onlyCoding=<boolean> --pathogenic=<boolean>]
\n options:\n --maffile=<file> maffile that has been annotated using vep & oncokb.
\n --outputname=<string> output string
\n --ACMG=<file> Table of ACMG variants
\n --pathwayList=<string> a List to return variants of a specific pathway (believed to be related to response) 
\n --pathwayFile=<file> a Table of pathway terms to search
\n --caddscore=<float> [default: 0] score for cadd cut-off
\n --gnomadcutoff=<float> [default: 0.1] 
\n --onlyCoding=<boolean> [default: TRUE]
\n --pathogenic=<boolean> [default: TRUE] filter to keep only variants which have a high CADD score OR non-benign ClinVar OR damaging POLYPHEN" -> doc

library("docopt", quietly = T)
opts <- docopt(doc)
print(opts)

suppressMessages(library(data.table, quietly = T))
suppressMessages(library(dplyr, quietly = T))

InputData=fread(opts$maffile, sep="\t")
print('RUN VUS SEARCH')
sprintf('Initial Number of variants %s', nrow(InputData))
## Remove the variants associated with 
rmidx1=which(InputData$CancerGeneCensus.Tier%in%c("Hallmark 1", "1", "Hallmark 2", "2")| InputData$CancerMutationCensus.Tier%in%c(1:3))
sprintf('Remove Cosmic related variants N=%s', length(rmidx1))
ex2=which(InputData$UserGeneList==1)
PWtable=read.csv(opts$pathwayFile)
nx=which(colnames(PWtable)==opts$pathwayList)
nx=setdiff(as.vector(PWtable[ ,nx]), "")
ex1=unique(unlist(sapply(nx, function(x) grep(x, InputData$HallmarkPathways))))
rmidx2=unique(c(ex1, ex2))
sprintf('Remove User gene list variants N=%s', length(rmidx2))
rmidx3=grep("drug_response|protective", InputData$ClinVar.Sig)
sprintf('Remove Clinvar protective variants N=%s', length(rmidx3))

ACMGtable=read.csv(opts$ACMG)
rmidx4=which(InputData$SYMBOL%in%ACMGtable$Gene)
sprintf('Remove all genes within ACMG list N=%s', length(rmidx4))

rmidx=unique(c(rmidx1, rmidx2, rmidx3))
idx=setdiff(1:nrow(InputData), rmidx)

Keep3=InputData[idx, ]

sprintf('Finding other VUS. Filter by population AF %s', as.numeric(opts$gnomadcutoff))
ax1=which(Keep3$AF_max<=as.numeric(opts$gnomadcutoff))
Keep3=Keep3[ax1,  ]
sprintf('%s variants kept after gnomad', nrow(Keep3))

if (opts$onlyCoding){
  print('Find variants only in coding areas')
  ax1=which(Keep3$ConsB==1)
  Keep3=Keep3[ax1,  ]
  sprintf('%s variants kept after filtering to coding regions', nrow(Keep3))
}

if (as.numeric(opts$caddscore)>0){
  sprintf('Find variants with cadd cutoff of %s', as.numeric(opts$caddscore))
  ax1=which(Keep3$CADD>=as.numeric(opts$caddscore))
  Keep3=Keep3[ax1,  ]
  sprintf('%s variants kept after cadd', nrow(Keep3))
}

#InputData$ProteinDomain!=" "

if (opts$pathogenic){
  print('filter based on predicted pathogenicity')
  rmT=which(Keep3$Pathogenicity==1)
  Keep3=Keep3[rmT, ]
  sprintf("%s variants after clinvar PolyPhen or CADD",nrow(Keep3))
}

Keep3=Keep3%>%select(!c("AF_max", "ConsB"))

sprintf("%s variants are VUS ", nrow(Keep3))

write.table(Keep3, file=paste(opts$outputname, "pathogenicVUS.filt.maf", sep=""), sep = "\t", row.names = F,  quote = F)


