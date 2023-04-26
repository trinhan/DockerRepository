#!/usr/bin/Rscript
"
This script takes an AnnotSV annotated file, annotated using the SummarizeAnnotSV.R script. The outputs are:
 * a full transcript file (filtered for the columns of interest, and by minimum VAF or CN)
 * a split transcript file (filtered for genes of interest - either from Cosmic or Pathway related, and by minimum VAF or CN)
 * acmg transcript file (filtered according to ACMG)
 * summary table - number of different variants (after filtering)

usage: \n FilterSVs.R [--tsv=<file> --outputname=<string> --mode=<string> --VAF=<float> --ACMGcutoff=<int> --ColsIDs=<file> --CNlow=<int> --CNhigh=<int>] 
\n options:
\n --tsv=<file> tsv from AnnotSV or funcotator
\n --outputname=<string> output string
\n --mode=<string> CNV or SV mode
\n --VAF=<float> [default: 0.25] VAF threshold to use
\n --ACMGcutoff=<int> [default: 4] for germline mode 
\n --CNlow=<int> [default: 1] is this in germline mode
\n --CNhigh=<int> [default: 3] is this in germline mode
\n --ColsIDs=<file> [default: NULL] Matrix of required column inputs and how to rename them " -> doc

library("docopt", quietly = T)
opts <- docopt(doc, help=TRUE, version='1.0.0')

  ############################################
  #1. Load in all required libraries silently
  ############################################
  
  suppressMessages(library(data.table, quietly = T))
  suppressMessages(library(matrixStats, quietly = T))
  
  ############################################
  #2. Read in data file, set the CN thresholds and separate into split and full
  ############################################
  message('Run FilterSV')
  AnnotSV=read.delim(opts$tsv)
  
  if (opts$mode=="SV"){
  AnnotSV=AnnotSV[which(AnnotSV$MeanVAF>=opts$VAF), ]  ## filter criteria?
  sprintf('There are %s rows in data after %s VAF filter', nrow(AnnotSV), opts$VAF)
  }else if (opts$mode=="CNV"){
  AnnotSV=AnnotSV[which(AnnotSV$CN<=opts$CNlow |AnnotSV$CN>=opts$CNhigh), ]  
  sprintf('There are %s rows in data with CN less than or equal to %s and higher than or equal to %s ', nrow(AnnotSV), opts$CNlow, opts$CNhigh)
  }
  AnnotSVFull=AnnotSV[AnnotSV$Annotation_mode=="full", ]
  AnnotSVSplit=AnnotSV[AnnotSV$Annotation_mode=="split", ]
  
  # stop catch if there is no data in AnnotSV
  if(nrow(AnnotSV)==0){
    message("There are not enough rows in the input SV file...")
  }
  
  ###########################################
  # 3. Create a summary table
  ###########################################
  message('creating summary table')
  gain=which(AnnotSVFull$SV_type=="DUP")
  homloss=which(AnnotSVFull$SV_type=="DEL")
  ins=which(AnnotSVFull$SV_type=="INS")
  inv=which(AnnotSVFull$SV_type=="INV")
  bnd=which(AnnotSVFull$SV_type=="BND")
  
  VariantSummary=c("Number of DUP/GAIN", "Number of DEL",
                   "Number of INS", "Number of INV",
                   "Number of BND")
  
  Nx=data.frame(Ngain=length(gain), Ndel=length(homloss), Nins=length(ins), Ninv=length(inv), Nbnd=length(bnd))
  SummTable=cbind(VariantSummary, N=t(Nx))

  ############################################
  #2B. Prepare the dictionary for renaming data
  ############################################
  if (opts$ColsID=="NULL"){
    SIDCol=read.csv("./annotFiles/SV_column_IDs.csv")
  } else {
    SIDCol=read.csv(opts$ColsID)
  }
    
  ############################################
  #3. Filter full length
  ############################################

  # Filter based on whether there is a GeneofInterest or Cosmic related gene
  # Add a tag to indicate the genes of interest
  if (opts$mode=="SV"){
    scolidx=which(SIDCol$Table%in%c("ALL", "FULL", "SV"))
  }else if (opts$mode=="CNV"){
    scolidx=which(SIDCol$Table%in%c("ALL", "FULL", "CNV"))
  }
  m1=match(SIDCol$Variable[scolidx], colnames(AnnotSVFull))
  S1=AnnotSVFull[, na.omit(m1)]
  colnames(S1)=SIDCol$ReName[scolidx[which(!is.na(m1))]]
  ## Reannotate based on chromosome
  S1$CHROM=factor(S1$CHROM, levels=c(1:22, "X", "Y"))
  S1=S1[order(S1$CHROM,S1$START), ]
  sprintf('Full length SVs associated with pathway/Cosmic: %s entries', nrow(S1))
  
  ############################################
  # 4. Filter split length
  ############################################
  #find the variables of interest and rename
  if (opts$mode=="SV"){
    scolidx=which(SIDCol$Table%in%c("ALL", "SPLIT", "SV"))
  }else if (opts$mode=="CNV"){
    scolidx=which(SIDCol$Table%in%c("ALL", "SPLIT", "CNV"))
  } 
  m1=match(SIDCol$Variable[scolidx], colnames(AnnotSVSplit))
  S2=AnnotSVSplit[, na.omit(m1)]
  colnames(S2)=SIDCol$ReName[scolidx[which(!is.na(m1))]]
  #########################################
  # !!modify this to include all genes
  ########################################
  idx=which(AnnotSVSplit$GenesOfInterest!=""|AnnotSVSplit$Cosmic!="")
  #try(if(length(idx)==0) warning("No genes of interest found"))
  #S2=S2[idx, ]
  #idx=length(unique(S2$GENE))
  
  #if(length(idx)>0){
  # also indicate whether to use a gene for a plot based on how long the CV is, or whether it crosses multiple intron/exon boundaries
    SpanReg=strsplit(S2$Location,"-")
    StartIdx=sapply(SpanReg, function(x) x[1])
    EndIdx=sapply(SpanReg, function(x) x[2])
    MultiInt=ifelse(StartIdx!=EndIdx, T, F)
    MultiInt[grep("exon", StartIdx)]=T
    MultiInt[grep("exon", EndIdx)]=T
    S2$MultiSpan=MultiInt
    S2$CHROM=factor(S2$CHROM, levels=c(1:22, "X", "Y"))
    S2=S2[order(S2$CHROM,S2$START), ]
  #}
    
  sprintf('Split length SVs total %s , and associated with pathway/Cosmic: %s entries', nrow(S2), length(idx))
  
  ###############################################
  # 5. Filter the ACMG genes
  ###############################################
  message('find the acmg genes')
  filt2=which(AnnotSVFull$ACMG_class>opts$ACMGcutoff)
  S3=S1[filt2, ]
  #colnames(S3)[c(1:3,7:8)]=c("CHROM", "START", "END", "Pop.Freq","PathogenicSource")
  sprintf('Filter acmg SVs with %s cutoff: %s entries', opts$ACMGcutoff,nrow(S3))
  
  ################################################
  # 6. Write to file, or return the list
  ###############################################
  message('write files to output')
  write.table(SummTable, file=paste(opts$outputname,".", opts$mode, ".SummaryTable.txt", sep=""), sep = "\t", row.names = F,  quote = F)
  write.table(S1, file=paste(opts$outputname,".", opts$mode, ".full.filt.maf", sep=""), sep = "\t", row.names = F,  quote = F)
  write.table(S2, file=paste(opts$outputname,".", opts$mode, ".split.filt.maf", sep=""), sep = "\t", row.names = F,  quote = F)
  write.table(S3, file=paste(opts$outputname,".", opts$mode, ".acmg.filt.maf", sep=""), sep = "\t", row.names = F,  quote = F)

