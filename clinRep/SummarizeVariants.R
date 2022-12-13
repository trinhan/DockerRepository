#!/usr/bin/Rscript
"usage: \n SummarizeVariants.R [--maffile=<file> --outputname=<string> --germline=<boolean> --caddscore=<float> --gnomadcutoff=<float> --Ncallerthresh=<int> --AddList=<file> --columnEntries=<file>]
\n options:\n --maffile=<file> maffile that has been annotated using vep & oncokb.
\n --outputname=<string> output string
\n --germline=<boolean> [default: TRUE] is this in germline mode
\n --caddscore=<float> [default: 20] score for cadd cut-off
\n --gnomadcutoff=<float> [default: 0.1]
\n --Ncallerthresh=<int> [default: 0] minimum number of callers required to output the variant
\n --AddList=<file> Additional user defined genes [default=NULL]
\n --columnEntries=<file> Import a table of columns to export and renamed output [default=NULL] " -> doc

library("docopt", quietly = T)
opts <- docopt(doc)
#SummarizeVariants(opts$maffile, opts$outputname, opts$germline, opts$caddscore, opts$gnomadcutoff, opts$exomesize, opts$pathwayList)
#SummarizeVariants=function(maffile, outputname, germline=T, caddscore=20, gnomadcutoff=0.1, exomesize=45, pathwayList,ACMG, AddList=NULL){
suppressMessages(library(data.table, quietly = T))
suppressMessages(library(dplyr, quietly = T))
  writeLines('Running SummarizeVariants \n Read in original File.. \n')
  mydnames=fread(opts$maffile, sep="\t", nrows=0)
  gnames=colnames(mydnames)[grep("nTot", colnames(mydnames))]
  gnames2=colnames(mydnames)[grep("^GT", colnames(mydnames))]
  gnames3=colnames(mydnames)[grep("VAF", colnames(mydnames))]
  gnames4=colnames(mydnames)[grep("nAlt", colnames(mydnames))]
  
  ColNames=read.csv(opts$columnEntries, header=F)
  
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
  
  newNames
  
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
    Test1$gnomAD.AF[which(is.na(Test1$gnomAD.AF))]=0
    Test1$MGRB.AF[which(is.na(Test1$MGRB.AF))]=0
    Test1$AF_max=ifelse(Test1$gnomAD.AF>Test1$MGRB.AF, Test1$gnomAD.AF, Test1$MGRB.AF)
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
  Nx=data.matrix(Nx)
  Ncall=rowSums(sign(Nx), na.rm = T)
  Ncallers=sapply(1:nrow(Nx2), function(x) paste(na.omit(Nx2[x, ]),collapse=","))
  head(Ncallers)
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
  
  print('Modify disease frequency summary')
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
  
  ## Increased cancer disease freq
  sprintf('Keep only variants called by more than %s callers', opts$Ncallerthresh)
  Test1=Test1[which(Test1$Ncall>opts$Ncallerthresh), ]

  
  ## Test whether in additional annot List
  Test1$UserGeneList=NA
  if (!is.null(opts$AddList)){
    ImmTable=read.csv(opts$AddList, header=F)
    ImmTable=ImmTable[ ,1]
    int2=which(Test1$SYMBOL%in%ImmTable)
    Test1$UserGeneList[int2]=1
  }
  
  ## Summary of variants which fits the cadd cutoffs, gnomad AF cutoffs and is in a coding region
  ConsequenceVals=c("missense", "nonsense", "frameshift", "splice", "UTR", "inframe", "stop", "NMD")
  Nxgrep=sapply(ConsequenceVals, function(x) grep(x, Test1$Consequence))
  Test1$ConsB=NA
  Test1$ConsB[unlist(Nxgrep)]=1
  
  # Add a potential filter based on pathogenicity
  Test1$Pathogenicity=NA
  # select based on clinVar
  select3=unique(unlist(sapply(c("pathogenic", "risk_factor", "drug_response", "protective"), function(x) grep(x, Test1$ClinVar.Sig, ignore.case=T))))
  length(select3)
  # select based on polyphen
  select4=grep("damaging", Test1$PolyPhen)
  length(select4)
  # select based on CADD
  if (as.numeric(opts$caddscore)>0){
    select5=which(Test1$CADD> as.numeric(opts$caddscore))
    length(select5)
    sprintf("%s variants after CADD filter" ,length(select5))
    idx=unique(c(select3, select4, select5))
  }else{
    idx=unique(c(select3, select4))
  }
  
  Test1$Pathogenicity[idx]=1
  
  write.table(Test1, file=paste(opts$outputname, "variantsAll.maf", sep=""), sep = "\t", row.names = F,  quote = F)
    
  # print('Finding ACMG variants')
  # ## ACMG variants
  # ACMGtable=read.csv(opts$ACMG)
  # Lx1=which(Test1$SYMBOL%in%ACMGtable$Gene)
  # strSplitA=strsplit(Test1$ClinVar.Sig, "&")
  # Nx2=sapply(strSplitA, function(x) length(which(c("pathogenic", "Pathogenic", "Likely_pathogenic", "likely_pathogenic")%in%x)))
  # Lx2=which(Nx2>0)
  # Lx3=intersect(Lx1, Lx2)
  # sprintf("%s variants in ACMG list", length(Lx3))
  # Keep0=Test1[Lx3, ]%>%select(!"AF_max")

  #if (length(Lx3)>0){
  #  Keep0=Test1[Lx3, ]
  #}else{
  #  Keep0=""
  #}
  
#   print('Finding cancer variants')
#   ## Known Variants associated with cancer?
#   bx1=which((Test1$CancerGeneCensus.Tier%in%c("Hallmark 1", "1", "Hallmark 2", "2") & Test1$CancerMutationCensus.Tier%in%c(1:3)))
#  # length(bx1)
#  # bx2=which(Test1$CMC.Cancer_Gene_Tier %in% c("Hallmark 1", "1", "Hallmark 2", "2"))
# #  table(Test1$CMC.Cancer_Gene_Tier)
#   #colnames(Test1)
#   sprintf('%s in Cosmic', length(bx1))
#   cx1=grep("Oncogenic", Test1$Oncokb.ONCOGENIC)
#   sprintf('%s in Oncokb', length(cx1))
#   mx1=unique(sort(c(bx1, cx1)))
#   sprintf('Total %s cancer variants', length(mx1))
#   Keep1=Test1[mx1, ] %>% select(!"AF_max")
#   
  # print('Finding drug assoc variants')
  # ## Known Variants associated with drug response or 
  # nx1=grep("drug_response|protective", Test1$ClinVar.Sig)
  # sprintf("%s variants realted to drug response or protective effect", length(nx1))
  # Keep2=Test1[nx1, ]%>%select(!"AF_max")
  
 

  # sprintf('Find genes involved in %s pathway', opts$pathwayList)
  # ## Pathways of Interest /opt/PathwayList.csv
  # PWtable=read.csv("/annotFiles/PathwayList.csv")
  # nx=which(colnames(PWtable)==opts$pathwayList)
  # nx=setdiff(as.vector(PWtable[ ,nx]), "")
  # ex1=unique(unlist(sapply(nx, function(x) grep(x, Test1$HallmarkPathways))))
  # sprintf( "%s variants pathway of interest" ,length(ex1))  
  # 
  # if (!is.null(opts$AddList)){
  #   ex2=which(Test1$UserGeneList==1)
  #   ex1=unique(c(ex1, ex2))
  #   sprintf( "%s variants in user defined list" ,length(ex2))
  #   }
  # 
  # Keep4=Test1[ex1, ]%>%select(!c("AF_max", "ConsB"))
  # 
  # # Refine this to include CADD high score if it does not affect a coding region
  # #lx1=which(is.na(Keep4$HGVSp))
  # if (opts$caddscore>0){
  #   print('keep variants with high cadd cutoff')
  #   select2=which(Keep4$CADD>as.numeric(opts$caddscore))
  # } else{
  #   print('keep variants with high cadd cutoff')
  #   select2=which(Keep4$ConsB==1)
  #  }
  # select3=unique(unlist(sapply(c("pathogenic", "risk_factor", "drug_response"), function(x) grep(x, Keep4$ClinSig))))
  # select4=grep("damaging", Keep4$PolyPhen)
  # rmT=unique(c(select2, select3, select4))
  # Keep4=Keep4[rmT, ]

  
  print('generate coding list based on gnomad, Cons and Pathogenicity')
  nidx=which(Test1$AF_max<=as.numeric(opts$gnomadcutoff) & Test1$ConsB==1 & Test1$Pathogenicity==1)
  sprintf("Summarize variants with gnomad %s Consequence coding, Expanded list contains %s variants", as.numeric(opts$gnomadcutoff), length(nidx))
  
  FilteredTable=Test1[nidx,  ]%>%select(!c("AF_max", "ConsB"))   
    ## Summary Stats
  DFValues=data.frame(Nvar=nrow(Test1))
  
  print('Finished! Writing out files')
  write.table(Test1, file=paste(opts$outputname, "variantsAll.maf", sep=""), sep = "\t", row.names = F,  quote = F)
  write.table(DFValues, file=paste(opts$outputname, "variantSummary.filt.maf", sep=""), sep = "\t", row.names = F,  quote = F)
  write.table(FilteredTable, file=paste(opts$outputname, "variantsCoding.filt.maf", sep=""), sep = "\t", row.names = F,  quote = F)
#}

