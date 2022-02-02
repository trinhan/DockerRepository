#!/usr/bin/Rscript
"usage: \n SummarizeVariants.R [--maffile=<file> --outputname=<string> --germline=<boolean> --caddscore=<integer> --gnomadcutoff=<float> --exomesize=<float> --pathwayList=<string> --ACMG=<file> --ImmuneList=<file>]
\n options:\n --maffile=<file> maffile that has been annotated using vep & oncokb.
\n --outputname=<string> output string
\n --germline=<boolean> [default: TRUE] is this in germline mode
\n --caddscore=<integer> [default: 20] score for cadd cut-off
\n --gnomadcutoff=<float> [default: 0.1]
\n --exomesize=<float> [default: 45] how big the indicated coding regions? For Mutational Burden calculations
\n --pathwayList=<string> a List to return variants of a specific pathway (believed to be related to response)
\n --ACMG=<file> Table of ACMG variants
\n --ImmuneList=<file> List of immune related genes [default=NULL]" -> doc

library("docopt")
opts <- docopt(doc)

SummarizeVariants=function(maffile, outputname, germline=T, caddscore=20, gnomadcutoff=0.1, exomesize=45, pathwayList,
                            ImmuneList=NULL){
  library(data.table)
  library(dplyr)
  mydnames=fread(maffile, sep="\t", nrows=0)
  gnames=colnames(mydnames)[grep("nTot", colnames(mydnames))]
  gnames2=colnames(mydnames)[grep("GT", colnames(mydnames))]
  gnames3=colnames(mydnames)[grep("VAF", colnames(mydnames))]
  InputData=fread(maffile, sep="\t", select=c("Hugo_Symbol", "CHROM","POS", "HGVSp_Short","HGVSc","EXON","Consequence","gnomAD_AF",
               "ONCOGENIC","pfam","pirsf", "CADD_PHRED", "ClinVar_CLNSIG","ClinVar_CLNDN", "ClinVar","CMC.CGC_TIER",
               "CMC.MUTATION_SIGNIFICANCE_TIER", "CMC.ONC_TSG", "CMC.DISEASE", "HallmarkPathways", gnames, gnames2, gnames3))
  
  
  rm1=c(grep("HLA", InputData$Hugo_Symbol))
  Test1=InputData[setdiff(1:nrow(InputData), rm1), 1:20 ]
  colnames(Test1)=c("SYMBOL", "CHROM","POS", "HGVSp","HGVSc","EXON","Consequence","gnomAD.AF",
                    "Oncokb.ONCOGENIC","pfam","pirsf", "CADD", "ClinVar.Sig","ClinVar.Disease", "ClinVar","CancerGeneCensus.Tier",
                    "CancerMutationCensus.Tier", "CMC.Oncogene", "Cancer.Disease.Freq", "HallmarkPathways")
  
  ## Rename the Test1
  Test1$ClinVar.Disease=gsub("not_specified", "", Test1$ClinVar.Disease)
  Test1$ClinVar.Disease=gsub("not_provided", "", Test1$ClinVar.Disease)
  Test1$ClinVar.Disease=gsub("&", " ", Test1$ClinVar.Disease)
  Test1$ProteinDomain=ifelse(is.na(Test1$pfam), Test1$pirsf, Test1$pfam)
  ##Test1=data.frame(Test1)
 
  Nx=InputData %>% select(gnames)
  Nx=data.frame(Nx)
  AvDepth=ceiling(rowMeans(Nx, na.rm=T))
  HeaderV=strsplit(colnames(Nx),"\\.")
  HeaderV=sapply(HeaderV, function(x) x[3])
  Nx2=sapply(1:ncol(Nx), function(x) ifelse(!is.na(Nx[,x ]), HeaderV[x], NA))
  Ncall=rowSums(sign(Nx), na.rm = T)
  Ncallers=sapply(1:nrow(Nx2), function(x) paste(na.omit(Nx2[x, ]),collapse=","))
  Test1$Ncallers=Ncallers[setdiff(1:nrow(InputData), rm1)]
  Test1$Ncall=Ncall[setdiff(1:nrow(InputData), rm1)]
  Test1$Depth=AvDepth[setdiff(1:nrow(InputData), rm1)]
  ## Do the same thing with the GT
  GT=InputData %>% select(gnames2)
  GT=data.frame(GT)
  m1=which(GT[ ,1]!=GT[,2])
  VAF=InputData %>% select(gnames3)
  VAF2=rowMeans(VAF, na.rm = T)
  
  Genotype=ifelse(!is.na(GT[ ,3]), GT[ ,3], ifelse(!is.na(GT[,2]), GT[,2], GT[,1]))
  Test1$VAF=VAF2[setdiff(1:nrow(InputData), rm1)]
  Test1$Genotype=Genotype[setdiff(1:nrow(InputData), rm1)]
  
  ## delete this later
  ###est1$Genotype=ifelse(Test1$VAF>0.85, "1/1", "0/1")
  
  ## Increased cancer disease freq
  Test1=Test1[which(Test1$Ncall>1), ]
  Test1$Cancer.Disease.Freq=gsub("\\/NS", "", Test1$Cancer.Disease.Freq)
  Test1$Cancer.Disease.Freq=gsub("\\/carcinoma\\/", " ", Test1$Cancer.Disease.Freq)
  Test1$Cancer.Disease.Freq=gsub("squamous_cell_carcinoma", "SCC", Test1$Cancer.Disease.Freq)
  Test1$Cancer.Disease.Freq=gsub("non_small_cell_carcinoma", "NSCC", Test1$Cancer.Disease.Freq)
  Test1$Cancer.Disease.Freq=gsub("=[0-9]*\\/[0-9]*", "", Test1$Cancer.Disease.Freq)
  Test1$Cancer.Disease.Freq=gsub("\\/", " ", Test1$Cancer.Disease.Freq)
  
  ## Test whether in additional annot List
  
  if (!is.null(ImmTable)){
    ImmTable=read.csv(ImmuneList)
    int2=which(Test1$SYMBOL%in%ImmTable$V1)
    Test1$ImmuneList=NA
    Test1$ImmuneList[int2]=1
  }
  
  
  print('finding cancer assoc variants')
  ## ACMG variants

  ACMGtable=read.csv(ACMGtable)
  Lx1=which(Test1$SYMBOL%in%ACMGtable$Gene)
  Lx2=grep("Pathogenic|Likely_pathogenic", Test1$ClinVar.Sig)
  Lx3=intersect(Lx1, Lx2)
  
  if (length(Lx3)>0){
    Keep0=Test1[Lx3, ]
  }else{
    Keep0=""
  }
  
  ## Known Variants associated with cancer?
  bx1=which(Test1$CancerGeneCensus.Tier%in%c("1", "2", "Hallmark") & Test1$CancerMutationCensus.Tier%in%c(1:3))
  cx1=grep("Oncogenic", Test1$Oncokb.ONCOGENIC)
  mx1=unique(sort(c(bx1, cx1)))
  Keep1=Test1[mx1, ]
  print('finding drug assoc variants')
  ## Known Variants associated with drug response or 
  nx1=grep("drug_response|protective", Test1$ClinVar.Sig)
  Keep2=Test1[nx1, ]
  print('finding other VUS')
  print(caddscore)
  ax1=which(Test1$CADD>caddscore & (Test1$gnomAD.AF<gnomadcutoff | is.na(Test1$gnomAD.AF)) & Test1$ProteinDomain!=" ")
  Keep3=Test1[setdiff(ax1, c(nx1, mx1)), ]
  print(nrow(Keep3))
  ## Pathways of Interest
  PWtable=read.csv("/opt/PathwayList.csv")
  nx=which(colnames(PWtable)==pathwayList)
  nx=setdiff(as.vector(PWtable[ ,nx]), "")
  ex1=unique(unlist(sapply(nx, function(x) grep(x, Test1$HallmarkPathways))))
  
  if (!is.null(ImmTable)){
    ex2=which(Test1$ImmuneList==1)
    ex1=unique(c(ex1, ex2))
    }
  
  Keep4=Test1[ex1, ]

  # Refine this to include CADD high score if it does not affect a coding region
  lx1=which(is.na(Keep4$HGVSp))
  rmT=which(Keep4$ClinVar.Sig!="Benign" | Keep4$CADD>caddscore)
  Keep4=Keep4[setdiff(rmT, lx1), ]
  
  ## Summary Stats
  DFValues=data.frame(Nvar=nrow(Test1), Cancer=length(mx1), DrugResp=length(nx1), Hallmark=nrow(Keep4),OtherVUS=nrow(Keep3))
  
  print('Write out files')
  write.table(Keep0, file=paste(outputname, "ACMG.filt.maf", sep=""), sep="\t", row.names = F, quote = F)
  write.table(Keep1, file=paste(outputname, "cancerVariants.filt.maf", sep=""), sep="\t", row.names = F, quote = F)
  write.table(Keep2, file=paste(outputname, "drugprotective.filt.maf", sep=""), sep="\t", row.names = F, quote = F)
  write.table(Keep3, file=paste(outputname, "pathogenicVUS.filt.maf", sep=""), sep="\t", row.names = F, quote = F)
  write.table(Keep4, file=paste(outputname, "HallmarkVUS.filt.maf", sep=""), sep="\t", row.names = F, quote = F)
  write.table(DFValues, file=paste(outputname, "variantSummary.filt.maf", sep=""), sep = "\t", row.names = F,  quote = F)
}

SummarizeVariants(opts$maffile, opts$outputname, opts$germline, opts$caddscore, opts$gnomadcutoff, opts$exomesize, opts$pathwayList)