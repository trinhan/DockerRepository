#!/usr/bin/Rscript
"usage: \n SummarizeVariants.R [--maffile=<file> --outputname=<string> --germline=<boolean> --caddscore=<integer> --gnomadcutoff=<float> --exomesize=<float> --Ncallerthresh=<int> --pathwayList=<string> --ACMG=<file> --AddList=<file>]
\n options:\n --maffile=<file> maffile that has been annotated using vep & oncokb.
\n --outputname=<string> output string
\n --germline=<boolean> [default: TRUE] is this in germline mode
\n --caddscore=<integer> [default: 20] score for cadd cut-off
\n --gnomadcutoff=<float> [default: 0.1]
\n --exomesize=<float> [default: 45] how big the indicated coding regions? For Mutational Burden calculations
\n --Ncallerthresh=<int> [default: 0] minimum number of callers required to output the variant
\n --pathwayList=<string> a List to return variants of a specific pathway (believed to be related to response)
\n --ACMG=<file> Table of ACMG variants
\n --AddList=<file> Additional user defined genes [default=NULL]" -> doc

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
  InputData=fread(opts$maffile, sep="\t", select=c("Hugo_Symbol","Existing_variation", "CHROM","POS", "HGVSp_Short","HGVSc","EXON","Consequence","gnomADg_AF",
               "mgrb_MGRB_AF","ONCOGENIC","pfam","pirsf", "CADD_PHRED","PolyPhen", "ClinVar_CLNSIG","ClinVar_CLNDN", "ClinVar","CMC.CGC_TIER", 
               "CMC.MUTATION_SIGNIFICANCE_TIER", "CMC.ONC_TSG", "CMC.DISEASE", "HallmarkPathways", gnames, gnames2, gnames3))
  
  
  rm1=c(grep("HLA", InputData$Hugo_Symbol))
  Test1=InputData[setdiff(1:nrow(InputData), rm1), 1:23 ]
  colnames(Test1)=c("SYMBOL", "rsID", "CHROM","POS", "HGVSp","HGVSc","EXON","Consequence","gnomAD.AF","MGRB.AF",
                    "Oncokb.ONCOGENIC","pfam","pirsf", "CADD", "PolyPhen", "ClinVar.Sig","ClinVar.Disease", "ClinVar","CancerGeneCensus.Tier",
                    "CancerMutationCensus.Tier", "CMC.Oncogene", "Cancer.Disease.Freq", "HallmarkPathways")
  
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
  Test1$VAF=round(VAF2[setdiff(1:nrow(InputData), rm1)], 2)
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
    ImmTable=read.csv(opts$AddList)
    int2=which(Test1$SYMBOL%in%ImmTable$V1)
    Test1$UserGeneList[int2]=1
  }
  

  print('Finding ACMG variants')
  ## ACMG variants
  ACMGtable=read.csv(opts$ACMG)
  Lx1=which(Test1$SYMBOL%in%ACMGtable$Gene)
  Lx2=grep("Pathogenic|Likely_pathogenic", Test1$ClinVar.Sig)
  Lx3=intersect(Lx1, Lx2)
  
  Keep0=Test1[Lx3, ]

  #if (length(Lx3)>0){
  #  Keep0=Test1[Lx3, ]
  #}else{
  #  Keep0=""
  #}
  
  print('Finding cancer variants')
  ## Known Variants associated with cancer?
  bx1=which(Test1$CancerGeneCensus.Tier%in%c("1", "2", "Hallmark") & Test1$CancerMutationCensus.Tier%in%c(1:3))
  cx1=grep("Oncogenic", Test1$Oncokb.ONCOGENIC)
  mx1=unique(sort(c(bx1, cx1)))
  Keep1=Test1[mx1, ]
  
  print('Finding drug assoc variants')
  ## Known Variants associated with drug response or 
  nx1=grep("drug_response|protective", Test1$ClinVar.Sig)
  Keep2=Test1[nx1, ]
  
  sprintf('Finding other VUS. Filter by gnomad AF %s and CADD score of %s', opts$gnomadcutoff, opts$caddscore )
  
  ax1=which(Test1$CADD>opts$caddscore & (Test1$gnomAD.AF<opts$gnomadcutoff | Test1$MGRB.AF<opts$gnomadcutoff | is.na(Test1$gnomAD.AF)) & Test1$ProteinDomain!=" ")
  Keep3=Test1[setdiff(ax1, c(nx1, mx1)), ]

  sprintf('Find genes involved in %s pathway', opts$pathwayList)
  ## Pathways of Interest
  PWtable=read.csv("/opt/PathwayList.csv")
  nx=which(colnames(PWtable)==opts$pathwayList)
  nx=setdiff(as.vector(PWtable[ ,nx]), "")
  ex1=unique(unlist(sapply(nx, function(x) grep(x, Test1$HallmarkPathways))))
  
  if (!is.null(opts$AddList)){
    ex2=which(Test1$UserGeneList==1)
    ex1=unique(c(ex1, ex2))
    }
  
  Keep4=Test1[ex1, ]

  # Refine this to include CADD high score if it does not affect a coding region
  lx1=which(is.na(Keep4$HGVSp))
  rmT=which(Keep4$ClinVar.Sig!="Benign" | Keep4$CADD>opts$caddscore)
  Keep4=Keep4[setdiff(rmT, lx1), ]
  
  ## Summary Stats
  DFValues=data.frame(Nvar=nrow(Test1), Cancer=length(mx1), DrugResp=length(nx1), Hallmark=nrow(Keep4),OtherVUS=nrow(Keep3))
  
  print('Finished! Writing out files')
  write.table(Keep0, file=paste(opts$outputname, "ACMG.filt.maf", sep=""), sep="\t", row.names = F, quote = F)
  write.table(Keep1, file=paste(opts$outputname, "cancerVariants.filt.maf", sep=""), sep="\t", row.names = F, quote = F)
  write.table(Keep2, file=paste(opts$outputname, "drugprotective.filt.maf", sep=""), sep="\t", row.names = F, quote = F)
  write.table(Keep3, file=paste(opts$outputname, "pathogenicVUS.filt.maf", sep=""), sep="\t", row.names = F, quote = F)
  write.table(Keep4, file=paste(opts$outputname, "PathwayVariants.filt.maf", sep=""), sep="\t", row.names = F, quote = F)
  write.table(DFValues, file=paste(opts$outputname, "variantSummary.filt.maf", sep=""), sep = "\t", row.names = F,  quote = F)
#}

