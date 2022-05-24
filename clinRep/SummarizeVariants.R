#!/usr/bin/Rscript
"usage: \n SummarizeVariants.R [--maffile=<file> --outputname=<string> --germline=<boolean> --caddscore=<float> --gnomadcutoff=<float> --exomesize=<float> --Ncallerthresh=<int> --pathwayList=<string> --ACMG=<file> --AddList=<file>]
\n options:\n --maffile=<file> maffile that has been annotated using vep & oncokb.
\n --outputname=<string> output string
\n --germline=<boolean> [default: TRUE] is this in germline mode
\n --caddscore=<float> [default: 20] score for cadd cut-off
\n --gnomadcutoff=<float> [default: 0.1]
\n --exomesize=<float> [default: 45] how big the indicated coding regions? For Mutational Burden calculations
\n --Ncallerthresh=<int> [default: 0] minimum number of callers required to output the variant
\n --pathwayList=<string> a List to return variants of a specific pathway (believed to be related to response)
\n --ACMG=<file> Table of ACMG variants
\n --AddList=<file> Additional user defined genes [default=NULL]
\n --columnEntries=<file> Import a table of columns to export and renamed output [default=NULL] " -> doc

library("docopt", quietly = T)
opts <- docopt(doc)
#SummarizeVariants(opts$maffile, opts$outputname, opts$germline, opts$caddscore, opts$gnomadcutoff, opts$exomesize, opts$pathwayList)
#SummarizeVariants=function(maffile, outputname, germline=T, caddscore=20, gnomadcutoff=0.1, exomesize=45, pathwayList,ACMG, AddList=NULL){
  library(data.table, quietly = T)
  library(dplyr, quietly = T)
  print('Read in original File')
  mydnames=fread(opts$maffile, sep="\t", nrows=0)
  gnames=colnames(mydnames)[grep("nTot", colnames(mydnames))]
  gnames2=colnames(mydnames)[grep("GT", colnames(mydnames))]
  gnames3=colnames(mydnames)[grep("VAF", colnames(mydnames))]
  gnames4=colnames(mydnames)[grep("nAlt", colnames(mydnames))]
  
  if (!is.null(opts$columnEntries)){
    print('column entries not null')
    ColNames=read.csv(opts$columnEntries, header=F)
  }else{
    print('column entries null')
    ColNames=read.csv("annotFiles/ColumnIDs.csv", header=F)
  }
  
  ## check which colnames mataches the original file
  midx=(match(ColNames$V1, colnames(mydnames)))
  mNames=colnames(mydnames)[na.omit(midx)]
  
  InputData=fread(opts$maffile, sep="\t", select=c(mNames, gnames, gnames2, gnames3, gnames4))
  #print(colnames(InputData))
  
  #print(ColNames$V2[which(!is.na(midx))])
  
  rm1=c(grep("HLA", InputData$Hugo_Symbol))
  Test1=InputData[setdiff(1:nrow(InputData), rm1), 1:length(mNames) ]
  newNames=ColNames$V2[which(!is.na(midx))]
  colnames(Test1)=ColNames$V2[which(!is.na(midx))]
  
  ## rename columns if they don't exist
  if ( !"ClinVar.Disease"%in%newNames & "ClinVar_Trait"%in%newNames){
    print('rename column clinvar')
    colnames(Test1)[which(newNames=="ClinVar_Trait")]="ClinVar.Disease"
  }
  
  if (!"ClinVar.Sig"%in%newNames & "ClinSig"%in%newNames){
     print('rename column ClinSig')
     colnames(Test1)[which(newNames=="ClinSig")]="ClinVar.Sig"
  }
  
  if (!"gnomAD.AF"%in%newNames & "gnomad_MAX_AF"%in%newNames){
    print('rename column gnomad')
    colnames(Test1)[which(newNames=="gnomad_MAX_AF")]="gnomAD.AF"
  }
  
  if (!"CADD"%in%newNames){
    print('reset CADD cut-off')
    opts$caddscore=0
  }
  
  if (!"MGRB.AF"%in%newNames & "gnomAD.AF"%in%colnames(Test1)){
    print('gnomad only')
    Test1$AF_max=Test1$gnomAD.AF
    Test1$AF_max[which(is.na(Test1$AF_max))]=0
  } else if ("MGRB.AF"%in%newNames & "gnomAD.AF"%in%colnames(Test1)){
    Test1$AF_max=ifelse(Test1$gnomAD.AF>Test1$MGRB.AF, Test1$gnomAD.AF, Test1$MGRB.AF)
    Test1$AF_max[which(is.na(Test1$AF_max))]=0
  }
  
  print('Parse ClinVar and Protein annotations')
  
  ## Rename the Test1
  Test1$ClinVar.Disease=gsub("not_specified", "", Test1$ClinVar.Disease)
  Test1$ClinVar.Disease=gsub("not_provided", "", Test1$ClinVar.Disease)
  Test1$ClinVar.Disease=gsub("&", " ", Test1$ClinVar.Disease)
  Test1$ProteinDomain=ifelse(is.na(Test1$pfam), Test1$pirsf, Test1$pfam)
  ##Test1=data.frame(Test1)
  print('Extract genotypes and callers')
  Nx=InputData %>% select(all_of(gnames))
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
  GT=InputData %>% select(all_of(gnames2))
  GT=data.frame(GT)
  ## collapse the genotypes if more than 1 column is found
  VAF=InputData %>% select(all_of(gnames3))
  VAF2=rowMeans(VAF, na.rm = T)
  Test1$VAF=round(VAF2[setdiff(1:nrow(InputData), rm1)], 3)
  ax=ncol(GT)
  if (ax>1){
  m1=which(GT[ ,1]!=GT[,2])
  Genotype=ifelse(!is.na(GT[ ,3]), GT[ ,3], ifelse(!is.na(GT[,2]), GT[,2], GT[,1]))
  } else{
  Genotype=GT[,1]  
  }
  Test1$Genotype=Genotype[setdiff(1:nrow(InputData), rm1)]
  ## delete this later
  ###est1$Genotype=ifelse(Test1$VAF>0.85, "1/1", "0/1")
  
  ## Increased cancer disease freq
  sprintf('Filter out variants called by more than %s callers', opts$Ncallerthresh)
  Test1=Test1[which(Test1$Ncall>opts$Ncallerthresh), ]
  
  Test1$Cancer.Disease.Freq=gsub("\\/NS", "", Test1$Cancer.Disease.Freq)
  Test1$Cancer.Disease.Freq=gsub("\\/carcinoma\\/", " ", Test1$Cancer.Disease.Freq)
  Test1$Cancer.Disease.Freq=gsub("squamous_cell_carcinoma", "SCC", Test1$Cancer.Disease.Freq)
  Test1$Cancer.Disease.Freq=gsub("non_small_cell_carcinoma", "NSCC", Test1$Cancer.Disease.Freq)
  Test1$Cancer.Disease.Freq=gsub("=[0-9]*\\/[0-9]*", "", Test1$Cancer.Disease.Freq)
  Test1$Cancer.Disease.Freq=gsub("\\/", " ", Test1$Cancer.Disease.Freq)
  
  ## Test whether in additional annot List
  Test1$UserGeneList=NA
  if (!is.null(opts$AddList)){
    ImmTable=read.csv(opts$AddList, header=F)
    ImmTable=ImmTable[ ,1]
    int2=which(Test1$SYMBOL%in%ImmTable)
    Test1$UserGeneList[int2]=1
  }
  
  print('Finding ACMG variants')
  ## ACMG variants
  ACMGtable=read.csv(opts$ACMG)
  Lx1=which(Test1$SYMBOL%in%ACMGtable$Gene)
  strSplitA=strsplit(Test1$ClinVar.Sig, "&")
  Nx2=sapply(strSplitA, function(x) length(which(c("pathogenic", "Pathogenic", "Likely_pathogenic", "likely_pathogenic")%in%x)))
  Lx2=which(Nx2>0)
  Lx3=intersect(Lx1, Lx2)
  sprintf("%s variants in ACMG list", length(Lx3))
  Keep0=Test1[Lx3, ]%>%select(!"AF_max")

  #if (length(Lx3)>0){
  #  Keep0=Test1[Lx3, ]
  #}else{
  #  Keep0=""
  #}
  
  print('Finding cancer variants')
  ## Known Variants associated with cancer?
  bx1=which(Test1$CancerGeneCensus.Tier%in%c("1", "2", "Hallmark") & Test1$CancerMutationCensus.Tier%in%c(1:3))
  sprintf('%s in Cosmic', length(bx1))
  cx1=grep("Oncogenic", Test1$Oncokb.ONCOGENIC)
  sprintf('%s in Oncokb', length(cx1))
  mx1=unique(sort(c(bx1, cx1)))
  sprintf('Total %s cancer variants', length(mx1))
  Keep1=Test1[mx1, ] %>% select(!"AF_max")
  
  print('Finding drug assoc variants')
  ## Known Variants associated with drug response or 
  nx1=grep("drug_response|protective", Test1$ClinVar.Sig)
  sprintf("%s variants realted to drug response or protective effect", length(nx1))
  Keep2=Test1[nx1, ]%>%select(!"AF_max")
  
  sprintf('Finding other VUS. Filter by gnomad AF %s and CADD score of %s', opts$gnomadcutoff, opts$caddscore )
  
  ConsequenceVals=c("missense", "nonsense", "frameshift", "splice", "UTR", "inframe")
  if (opts$caddscore>0){
  ax1=which(Test1$CADD>as.numeric(opts$caddscore) & Test1$AF_max<opts$gnomadcutoff & Test1$ProteinDomain!=" ")
  } else {
  Nxgrep=sapply(ConsequenceVals, function(x) grep(x, Test1$Consequence))
  Test1$ConsB=NA
  Test1$ConsB[unlist(Nxgrep)]=1
  ax1=which(Test1$AF_max<opts$gnomadcutoff & Test1$ProteinDomain!=" " & 
              Test1$ConsB==1)
  }
  sprintf("%s variants are VUS ", length(ax1))
  Keep3=Test1[setdiff(ax1, c(nx1, mx1)),  ]%>%select(!"AF_max")

  sprintf('Find genes involved in %s pathway', opts$pathwayList)
  ## Pathways of Interest /opt/PathwayList.csv
  PWtable=read.csv("annotFiles/PathwayList.csv")
  nx=which(colnames(PWtable)==opts$pathwayList)
  nx=setdiff(as.vector(PWtable[ ,nx]), "")
  ex1=unique(unlist(sapply(nx, function(x) grep(x, Test1$HallmarkPathways))))
  sprintf( "%s variants pathway of interest" ,length(ex1))  

  if (!is.null(opts$AddList)){
    ex2=which(Test1$UserGeneList==1)
    ex1=unique(c(ex1, ex2))
    sprintf( "%s variants in user defined list" ,length(ex2))
    }
  
  Keep4=Test1[ex1, ]%>%select(!"AF_max")

  # Refine this to include CADD high score if it does not affect a coding region
  lx1=which(is.na(Keep4$HGVSp))
  if (opts$caddscore>0){
    print('keep variants with high cadd cutoff')
    rmT=which(Keep4$ClinVar.Sig!="Benign" | Keep4$CADD>as.numeric(opts$caddscore))
  } else{
    print('keep variants with high cadd cutoff')
    rmT=which(Keep4$ClinVar.Sig!="Benign" & Keep4$ConsB==1)
  }
  Keep4=Keep4[setdiff(rmT, lx1),]
  
    
  ## Summary of variants which fits the cadd cutoffs, gnomad AF cutoffs and is in a coding region
  print('generate coding list')
  if (opts$caddscore>0) {
    nidx=which(Test1$AF_max<opts$gnomadcutoff & Test1$CADD>as.numeric(opts$caddscore) & Test1$ConsB==1)
    sprintf("Summarize variants withcadd %s, gnomad %s Consequence coding, Expanded list contains %s variants", opts$caddscore,opts$gnomadcutoff, length(nidx))
    
  }else{
    nidx=which(Test1$AF_max<opts$gnomadcutoff & Test1$ConsB==1)
    sprintf("Summarize variants with gnomad %s Consequence coding, Expanded list contains %s variants", opts$gnomadcutoff, length(nidx))
  }
  FilteredTable=Test1[nidx,  ]%>%select(!"AF_max")    
    ## Summary Stats
  DFValues=data.frame(Nvar=nrow(Test1), Cancer=length(mx1), DrugResp=length(nx1), Hallmark=nrow(Keep4), OtherVUS=nrow(Keep3), ACMG=nrow(Keep0), Filtered=length(nidx))
    
  
  print('Finished! Writing out files')
  write.table(Keep0, file=paste(opts$outputname, "ACMG.filt.maf", sep=""), sep="\t", row.names = F, quote = F)
  write.table(Keep1, file=paste(opts$outputname, "cancerVariants.filt.maf", sep=""), sep="\t", row.names = F, quote = F)
  write.table(Keep2, file=paste(opts$outputname, "drugprotective.filt.maf", sep=""), sep="\t", row.names = F, quote = F)
  write.table(Keep3, file=paste(opts$outputname, "pathogenicVUS.filt.maf", sep=""), sep="\t", row.names = F, quote = F)
  write.table(Keep4, file=paste(opts$outputname, "PathwayVariants.filt.maf", sep=""), sep="\t", row.names = F, quote = F)
  write.table(DFValues, file=paste(opts$outputname, "variantSummary.filt.maf", sep=""), sep = "\t", row.names = F,  quote = F)
  write.table(FilteredTable, file=paste(opts$outputname, "variantsCoding.filt.maf", sep=""), sep = "\t", row.names = F,  quote = F)
#}

