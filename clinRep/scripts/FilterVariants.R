#!/usr/bin/Rscript
"usage: 
FilterVariants.R ( --maffile=<file> --scoringRubrik=<file>)  [--outputname=<string> --ACMG=<file> --pathwayList=<string> --pathwayFile=<file> --caddscore=<float> --gnomadcutoff=<float> --onlyCoding=<boolean> --pathogenic=<boolean> ] 

\n options:
\n --maffile=<file> maffile that has been annotated using vep & oncokb.
\n --scoringRubrik=<file> Input file with weightings to score the inputs 
\n --outputname=<string> output string
\n --ACMG=<file> Table of ACMG variants
\n --pathwayList=<string> a List to return variants of a specific pathway (believed to be related to response) 
\n --pathwayFile=<file> a Table of pathway terms to search
\n --caddscore=<float> [default: 10] score for cadd cut-off
\n --gnomadcutoff=<float> [default: 0.1] 
\n --onlyCoding=<boolean> [default: TRUE]
\n --pathogenic=<boolean> [default: TRUE] filter to keep only variants which have a high CADD score OR non-benign ClinVar OR damaging POLYPHEN" -> doc

library("docopt", quietly = T)
opts <- docopt(doc, help=TRUE, version='1.0.0')
print(opts)

############################################
#1. Load in all required libraries silently
############################################
suppressMessages(library(data.table, quietly = T))
suppressMessages(library(dplyr, quietly = T))

#####################################################################
#2. Read in the maffile, and set-up the scoring matrix
#####################################################################  
InputData=fread(opts$maffile, sep="\t")
sprintf('Initial Number of variants %s', nrow(InputData))

ScorMat=matrix(0, nrow=nrow(InputData), ncol=12)
colnames(ScorMat)=c("gnomad", "coding", "cadd", "pathogenic", "clinvarpath", "clinvardrug", "ACMG", "Cosmic", "Oncokb", 
                    "Usergene","Hallmark", "Pathway")
ScorMat=data.frame(ScorMat)

#####################################################################
#1-4. Filter according to gnomad, only coding, cadd score, pathogenic
#####################################################################  

if (!is.null(opts$gnomad)){
  ax1=which(InputData$AF_max<=as.numeric(opts$gnomadcutoff))
  ScorMat$gnomad[ax1]=1
  sprintf('%s variants kept after gnomad', length(ax1))
}

if (opts$onlyCoding){
  message('Find variants only in coding areas')
  ax1=which(InputData$ConsB==1)
  ScorMat$coding[ax1]=1
  sprintf('%s variants kept after filtering to coding regions', length(ax1))
}

if (as.numeric(opts$caddscore)>0){
  sprintf('Find variants with cadd cutoff of %s', as.numeric(opts$caddscore))
  ax1=which(InputData$CADD>=as.numeric(opts$caddscore))
  ScorMat$cadd[ax1]=1
  sprintf('%s variants kept after cadd', length(ax1))
}

if (opts$pathogenic){
  ax1=which(InputData$Pathogenicity==1)
  ScorMat$pathogenic[ax1]=1
  sprintf("%s variants after clinvar PolyPhen or CADD", length(ax1))
}

#####################################################################
#5. Filter according ClinVar annotation
#####################################################################  

# filter based on clinvar pathogenic annotations

  strSplitA=strsplit(InputData$ClinVar.Sig, "&")
  Nx2=sapply(strSplitA, function(x) length(which(c("pathogenic", "Pathogenic", "Likely_pathogenic", "likely_pathogenic")%in%x)))
  Lx2=which(Nx2>0)
  sprintf("%s variants are pathogenic", length(Lx2))
  ScorMat$clinvarpath[Lx2]=1

# filter based on scoring annotations
  
  nx1=grep("drug_response|protective", InputData$ClinVar.Sig)
  sprintf("%s variants related to drug response or protective effect", length(nx1))
  ScorMat$clinvardrug[nx1]=1

#####################################################################
#6. Gene in ACMG list
#####################################################################

## ACMG variants
ACMGtable=read.csv(opts$ACMG)
Lx1=which(InputData$SYMBOL%in%ACMGtable$Gene)
sprintf('%s variants in ACMG list of genes', length(Lx1))
if (length(Lx1)>0){
  ScorMat$ACMG[Lx1]=1
}

#####################################################################
#7. Filter based on whether the variant is in Cosmic Mut or Gene List
#####################################################################  

  message('Finding cancer variants')
  ## Known Variants associated with cancer?
  # Find variants with Cancer Gene Census Tier in 1:2
  bx1=which(InputData$CancerGeneCensus.Tier%in%c("Hallmark 1", "1", "Hallmark 2", "2")) 
  # Find variants with Cancer Gene Census Mut Tier in 1:3
  cx1 = which(InputData$CancerMutationCensus.Tier%in%c(1:3))
  sprintf('%s in Cosmic Gene Census', length(bx1))
  sprintf('%s in Cosmic Mut Census', length(cx1))
  clist=unique(c(bx1, cx1))
  ScorMat$Cosmic[clist]=1
  
#####################################################################
#8. Filter based on presence in oncokb
#####################################################################  

dx1=grep("Oncogenic", InputData$Oncokb.ONCOGENIC)
sprintf('%s in Oncokb', length(dx1))
if (length(dx1)>0){
  ScorMat$Oncokb[dx1]=1
}
ScorMat$Cancer=ifelse(ScorMat$Oncokb==1|ScorMat$Cosmic==1, 1, 0)

#####################################################################
#9. Filter based on User Gene List and/or Pathway
#####################################################################  
# add user gene to matrix
ex2=which(InputData$UserGeneList==1)
if (length(ex2)>0){
  ScorMat$Usergene[ex2]=1
}
# add pathway to matrix
PWtable=read.csv(opts$pathwayFile)
nx=which(colnames(PWtable)==opts$pathwayList)
nx=setdiff(as.vector(PWtable[ ,nx]), "")

ex1=unique(unlist(sapply(nx, function(x) grep(x, InputData$HallmarkPathways))))
if (length(ex1)>0){
  ScorMat$Hallmark[ex1]=1
}
ScorMat$Pathway=ifelse(ScorMat$Usergene==1|ScorMat$Hallmark==1, 1, 0)


#####################################################################
#11. Upload the scoring matrix and Filter
#####################################################################  

# read in the weighting matrix
MatrixIn=read.csv(opts$scoringRubrik, row.names = 1)

for (i in 1:ncol(MatrixIn)){
  # for the column of interest, match the variable with column names
  t1=which(!is.na(MatrixIn[ ,i]))
  Nlist=rownames(MatrixIn)[t1]
  a1=as.matrix(ScorMat[ ,Nlist])
  b1=as.numeric(MatrixIn[t1 ,i])
  TempMat=a1%*%b1
  nx=length(which(b1>0))
  # Take the rowsum and find the max to identify the gene list
  Glist=which(TempMat==nx)
  Keep3=InputData[Glist, ]
  Keep3=Keep3%>%select(!c("AF_max", "ConsB"))
  message(sprintf("%s variants are %s ", nrow(Keep3), colnames(MatrixIn)[i]))
  write.table(Keep3, file=paste(opts$outputname, colnames(MatrixIn)[i],".filt.maf", sep=""), sep = "\t", row.names = F,  quote = F)
}

ScorMat2=cbind(InputData[ ,1:5], ScorMat)
write.table(ScorMat2, file=paste(opts$outputname, ".scoringMatrix.txt", sep=""), sep = "\t", row.names = F,  quote = F)

message('FILTER VARIANTS COMPLETE')




