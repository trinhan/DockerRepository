#!/usr/bin/Rscript
"usage: \n FindPathway.R [--maffile=<file> --outputname=<string> --pathwayFile=<file> --pathwayList=<string> --gnomadcutoff=<float> --onlyCoding=<string> --caddscore=<float>  --pathogenic=<string>]
\n options:\n --maffile=<file> maffile that has been annotated using vep & oncokb.
\n --outputname=<string> output string
\n --pathwayFile=<file> a Table of pathway terms to search
\n --pathwayList=<string> a List to return variants of a specific pathway (believed to be related to response) 
\n --gnomadcutoff=<float> [default: 0.1] 
\n --onlyCoding=<string> [default: true]
\n --caddscore=<float> [default: 0] score for cadd cut-off
\n --pathogenic=<string> [default: false] filter to keep only variants which have a high CADD score OR non-benign ClinVar OR damaging POLYPHEN" -> doc

library("docopt", quietly = T)
opts <- docopt(doc)
print(opts)

suppressMessages(library(data.table, quietly = T))
suppressMessages(library(dplyr, quietly = T))
print('RUN PATHWAY FILTER')
InputData=fread(opts$maffile, sep="\t")

sprintf('Find genes involved in %s pathway', opts$pathwayList)
print('Search based on pathway names')
## Pathways of Interest /opt/PathwayList.csv
PWtable=read.csv(opts$pathwayFile)
nx=which(colnames(PWtable)==opts$pathwayList)
nx=setdiff(as.vector(PWtable[ ,nx]), "")
ex1=unique(unlist(sapply(nx, function(x) grep(x, InputData$HallmarkPathways))))
sprintf( "%s variants pathway of interest" ,length(ex1))  
print('Search based on user gene list')
ex2=which(InputData$UserGeneList==1)
sprintf( "%s variants in user defined list" ,length(ex2))
Keep4=InputData[unique(c(ex1, ex2)), ]

sprintf('Filter based on gnomad cutoff of %s', as.numeric(opts$gnomadcutoff))
select1=which(Keep4$AF_max<=as.numeric(opts$gnomadcutoff))
Keep4=Keep4[select1, ]
sprintf("%s variants after gnomad filter" ,nrow(Keep4))

if (opts$onlyCoding){
  print('filter based on predicted pathogenicity')
  idx=which(Keep4$ConsB==1)
  Keep4=Keep4[idx, ]
  sprintf("%s variants in coding regions", nrow(Keep4))
}

# Refine this to include CADD high score if it does not affect a coding region
if (as.numeric(opts$caddscore>0)){
  print('filter based on cadd score')
  select2=which(Keep4$CADD>as.numeric(opts$caddscore))
  Keep4=Keep4[select2, ]
  sprintf("%s variants after CADD filter" ,nrow(Keep4))
} 

if (opts$pathogenic){
print('filter based on predicted pathogenicity')
rmT=which(Keep4$Pathogenicity==1)
Keep4=Keep4[rmT, ]
sprintf("%s variants after clinvar and PolyPhen",nrow(Keep4))
}

Keep4=Keep4%>%select(!c("AF_max"))
print('FINISHED PATHWAY SEARCH')
write.table(Keep4, file=paste(opts$outputname, "PathwayVariants.filt.maf", sep=""), sep="\t", row.names = F, quote = F)

