#!/usr/bin/Rscript
"usage: \n SummarizeFilterFunco.R [--tsv=<file> --outputname=<string> --MSigDB=<file> --GTex=<file> --CosmicList=<file> --AddList=<file> --pathwayTerm=<string> --pathwayList=<file> --Tissue=<string> --AFthreshold=<float> ]
\n options:
\n --tsv=<file> tsv from AnnotSV or funcotator
\n --outputname=<string> output string
\n --MSigDB=<file> File containing MsigDB gene set information
\n --GTex=<file> GTex Expression Data
\n --CosmicList=<file> File of Cosmic Related Genes
\n --AddList=<file> List of additional gene lists [default: NULL]
\n --pathwayTerm=<string>
\n --pathwayList=<file>
\n --Tissue=<string> Tissue type 
\n --AFthreshold=<float> [default: 0.1] Max pop frequency value for common variants" -> doc

library("docopt", quietly = T)
opts <- docopt(doc, help=TRUE, version='1.0.0')

############################################
#1. Load in all required libraries silently
############################################

suppressMessages(library(data.table, quietly = T))
suppressMessages(library(dplyr, quietly = T))
suppressMessages(library(GSEABase, quietly = T))
suppressMessages(library(matrixStats, quietly = T))

############################################
#2. Read in data file and Filter pass tags if available
############################################

InputData=read.delim(opts$tsv)
print(paste0("Number of rows:", nrow(InputData)))

############################################
#6. Extract the exact genes to search pathways - MSIGSB
############################################
print('Annotate implicated genes with pathway data')  

Genes=strsplit(InputData$genes, ",")
PathInH=getGmt(con=opts$MSigDB, geneIdType=SymbolIdentifier(),
               collectionType=BroadCollection(category="h"))
Mx2=geneIds(PathInH)
GS=unique(unlist(Genes))
GS2=paste("^", GS, "$", sep="")
Ids2=names(PathInH)
Ids2=gsub("HALLMARK_", "", Ids2)
names(Mx2)=Ids2
Mx3=stack(Mx2)
Tx2=sapply(GS, function(x) as.character(Mx3$ind[which(Mx3$values==x)]))
Nm2=sapply(Genes, function(x) match(x, names(Tx2)))
Nm3=sapply(Nm2, function(x) paste(unique(as.character(Tx2[x])), collapse=" "))
Nm3=gsub( "character\\(0\\)", "",Nm3)
Nm3=gsub("c\\(\"", "", Nm3)
Nm3=gsub("\", \"", " ", Nm3)
Nm3=gsub("\")", "", Nm3)
InputData$Pathways=Nm3

############################################
#7. Add user information on genes of interest and specific pathways
############################################

## also include information from the AddList
if (opts$AddList!="NULL"){
  print('Annotate with input gene list')
  GL=read.delim(opts$AddList, header=F)
  GL=GL[ ,1]
  Nx2=sapply(Genes, function(x) paste(x[which(x%in%GL)], collapse=", "))
  InputData$GenesOfInterest=Nx2
}else{
  InputData$GenesOfInterest=NA
}

# Based on user annotation
sprintf('Find genes involved in %s pathway', opts$pathwayTerm)
## Pathways of Interest
PWtable=read.csv(opts$pathwayList)
nx=which(colnames(PWtable)==opts$pathwayTerms)
nx=setdiff(as.vector(PWtable[ ,nx]), "")
ex1=unique(unlist(sapply(nx, function(x) grep(x, InputData$Pathways))))
InputData$GenesOfInterest[ex1]=ifelse(is.na(InputData$GenesOfInterest[ex1]), opts$pathwayTerm, paste(InputData$GenesOfInterest[ex1], opts$pathwayTerm))

############################################
#8.Add GTex data
############################################

sprintf('Annotate with GTex data, z-scores of %s relative to all other tissues', opts$Tissue)
# Also annotate these Genes according to GTex data
## "~/Downloads/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz"
GTex=read.delim(opts$GTex, skip=2)
l1=match(GS, GTex$Description)
GTex=GTex[na.omit(l1),  ]
rownames(GTex)=GS[which(!is.na(l1))]
gtissue=grep(opts$Tissue, colnames(GTex), ignore.case = T)
GZScore=(GTex[ ,gtissue]-rowMeans(GTex[ ,-c(1:2)]))/rowSds(data.matrix(GTex[ ,-c(1:2)]))
if (length(gtissue)>1){
  GZscore2=rowMaxs(data.matrix(GZScore), na.rm=T)
  names(GZscore2)=rownames(GZScore)
  GZScore=GZscore2
}
idx1=which(GZScore>0, arr.ind = T)
GTexNames=unique(names(idx1))
GTex2=GZScore[match(GTexNames, names(GZScore))]
GTex3=paste(GTexNames, "(", round(GTex2, 1), ")", sep="")

MatchCase=sapply(Genes, function(x) paste(na.omit(GTex3[match(x, GTexNames)]), collapse=" "))
## join all the GeneNames >> and add to the table of interest
InputData$GTex=MatchCase

############################################
#9.Add Cosmic Data
############################################

print('Add Cosmic Data')
CData=read.csv(opts$CosmicList)
Nx2=sapply(Genes, function(x) paste(x[which(x%in%CData$Gene.Symbol)], collapse=", "))
InputData$Cosmic=Nx2


############################################
#13. Print the raw table to file
############################################

print('Finished! Write raw formated output table')
print(paste0("Number of Annotated CNV rows:", nrow(InputData)))
write.table(InputData, file=paste(opts$outputname, ".CNV.formated.tsv", sep=""), sep="\t", row.names = F, quote = F)
