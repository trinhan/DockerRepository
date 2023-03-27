#!/usr/bin/Rscript
"
This script takes an AnnotSV annotated file, annotated using the SummarizeAnnotSV.R script. The outputs are:
 * a full transcript file (filtered for the columns of interest, and by minimum VAF)
 * a split transcript file (filtered for genes of interest - either from Cosmic or Pathway related, and by minimum VAF)
 * acmg transcript file (filtered according to ACMG)
 * summary table - number of different variants (after filtering)

usage: \n FilterSVs.R [--tsv=<file> --germline=<boolean> --VAF=<float> --ACMGcutoff=<int>]
\n options:
\n --tsv=<file> tsv from AnnotSV or funcotator
\n --outputname=<string> output string
\n --germline=<boolean> [default: TRUE] is this in germline mode
\n --VAF=<float> [default: 1] VAF threshold to use
\n --ACMGcutoff=<int> [default: 4] for germline mode " -> doc

library("docopt", quietly = T)
opts <- docopt(doc, help=TRUE, version='1.0.0')

FilterSVs=function(opts){
  ############################################
  #1. Load in all required libraries silently
  ############################################
  
  suppressMessages(library(data.table, quietly = T))
  suppressMessages(library(matrixStats, quietly = T))
  #load("./annotFiles/ensembl_granges_genes_290722_chr1_22_XYMT.RData")
  
  ############################################
  #2. Read in data file, set the CN thresholds and separate into split and full
  ############################################
  
  AnnotSV2=read.delim(opts$tsv)
  AnnotSV2=AnnotCNV[which(AnnotSV2$MeanVAF>=opts$VAF), ]  ## filter criteria?
  AnnotSVFull=AnnotCNV2[AnnotCNV2$Annotation_mode=="full", ]
  AnnotSVSplit=AnnotCNV2[AnnotCNV2$Annotation_mode=="split", ]
    
  
  ###########################################
  # 3. Create a summary table
  ###########################################
  gain=which(AnnotSVFull$SV_type=="DUP")
  homloss=which(AnnotSVFull$SV_type=="DEL")
  ins=which(AnnotSVFull$SV_type=="INS")
  inv=which(AnnotSVFull$SV_type=="INV")
  bnd=which(AnnotSVFull$SV_type=="BND")
  
  VariantSummary=c("Number of DUP/GAIN SVs ", "Number of DEL SVs ",
                   "Number of INS SVs ", "Number of INV SVs ",
                   "Number of BND SVs ")
  
  Nx=data.frame(Ngain=length(gain), Ndel=length(homloss), Nins=length(ins), Ninv=length(inv), Nbnd=length(bnd))
  SummTable= cbind(VariantSummary, t(Nx))
  
  ############################################
  #3. Filter full length
  ############################################
  Scolnames=c("Gene_name","AnnotSV_ID", "SV_type", "SV_length", "ACMG_class", "Location", "AF", "HI", "TS","OMIM_phenotype", "OMIM_morb",
              "RE_gene", "Dist_nearest_SS", "TAD_coordinate",
              "Bsource", "Psource", "GnomAD_pLI", "ExAC_pLI", "SRVAF","PRVAF",
              "MeanVAF", "Depth","GenesOfInterest","Cosmic", "Pathways","GTex", "SV_chrom","SV_start","SV_end"  )
  # Filter based on whether there is a GeneofInterest or Cosmic related gene
  # Add a tag to indicate the genes of interest
  S1=AnnotCNVFull[, Scolnames]
  colnames(S1)=c("GENE","SV_ID", "TYPE", "SV_length", "ACMG", "Location", "Pop.Freq", "Haploinsufficiency", "Triplosensitivity", "OMIM", "OMIM_morbid", "Reg.Element.Genes", 
                 "Nearest.Splice.Site.bp", "TAD", "BenignSource", "PathogenicSource", "gnomad.pLI", "exac.pLI", "SplitReadVAF", "PairedReadVAF", "MeanVAF", "Depth",
                 "GenesOfInterest", "Cosmic", "Pathways", "GTex", "CHROM" ,"START", "END")
  S1$CHROM=factor(S1$CHROM, levels=c(1:22, "X", "Y"))
  S1=S1[order(S1$CHROM,S1$START), ]
  
  ############################################
  # 4. Filter split length
  ############################################
  Scolnames=c("SV_chrom", "SV_start", "SV_end","Gene_name", "Pathways","CN", "GenesOfInterest","Cosmic","Pheno", "ACMG_class", "Location", "AF", "HI", "TS","OMIM_phenotype", "OMIM_morb",
              "RE_gene", "Dist_nearest_SS", "TAD_coordinate", "Bsource", "Psource", "GnomAD_pLI", "ExAC_pLI", "GTex")       
  # Filter based on whether there is a GeneofInterest or Cosmic related gene
  idx=which(AnnotCNVSplit$GenesOfInterest!=""|AnnotCNVSplit$Cosmic!="")
  S2=AnnotCNVSplit[idx, Scolnames]
  colnames(S2)=c("CHROM","START","END", "Gene_name","Pathways","CN","PathwayRelated","Cosmic", "Pheno","ACMG_class","Location","Pop.Freq","Haploinsufficiency","Triplosensitivity",
                 "OMIM","OMIM_morbid","Reg.Element.Genes","Nearest.Splice.Site.bp", "TAD","BenignSource","PathogenicSource","gnomad.pLI","exac.pLI","GTex")   
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
  
  ###############################################
  # 5. Filter the ACMG genes
  ###############################################
  Scolnames=c("SV_chrom", "SV_start", "SV_end", "CN", "Pheno", "ACMG_class", "AF", "Psource",
              "Gene_name","Cosmic", "Pathways","GTex","GenesOfInterest" )
  filt2=which(AnnotCNVFull$ACMG_class>opts$ACMGcutoff)
  S3=AnnotCNVFull[filt2, Scolnames]
  colnames(S3)[c(1:3,7:8)]=c("CHROM", "START", "END", "Pop.Freq","PathogenicSource")
  
  ################################################
  # 6. Write to file, or return the list
  ###############################################
  return(list(SummTable=SummTable,full=S1, split=S2, acmg=S3))
  
}