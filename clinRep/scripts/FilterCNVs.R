#!/usr/bin/Rscript
"usage: \n FilterCNVs.R [--tsv=<file> --germline=<boolean> --CNlow=<int> --CNhigh=<int> --ACMGcutoff=<int>]
\n options:
\n --tsv=<file> tsv from AnnotSV or funcotator
\n --outputname=<string> output string
\n --germline=<boolean> [default: TRUE] is this in germline mode
\n --CNlow=<int> [default: 1] is this in germline mode
\n --CNhigh=<int> [default: 3] is this in germline mode
\n --ACMGcutoff=<int> [default: 4] for germline mode " -> doc

library("docopt", quietly = T)
opts <- docopt(doc, help=TRUE, version='1.0.0')

FilterCNVs=function(opts){
############################################
#1. Load in all required libraries silently
############################################

suppressMessages(library(data.table, quietly = T))
suppressMessages(library(matrixStats, quietly = T))
#load("./annotFiles/ensembl_granges_genes_290722_chr1_22_XYMT.RData")

############################################
#2. Read in data file, set the CN thresholds and separate into split and full
############################################

AnnotCNV=read.delim(opts$tsv)
AnnotCNV2=AnnotCNV[which(AnnotCNV$CN<opts$CNlow |AnnotCNV$CN>opts$CNhigh), ]
AnnotCNVFull=AnnotCNV2[AnnotCNV2$Annotation_mode=="full", ]
AnnotCNVSplit=AnnotCNV2[AnnotCNV2$Annotation_mode=="split", ]

############################################
#3. Filter full length
############################################
Scolnames=c("SV_chrom", "SV_start", "SV_end", "CN", "Pheno", "ACMG_class", "Location", "AF", "HI", "TS","OMIM_phenotype", "OMIM_morb",
            "RE_gene", "Dist_nearest_SS", "TAD_coordinate",
            "Bsource", "Psource", "GnomAD_pLI", "ExAC_pLI","GenesOfInterest","Cosmic", "Pathways","GTex", "Gene_name")
# Filter based on whether there is a GeneofInterest or Cosmic related gene
# Add a tag to indicate the genes of interest
T1=AnnotCNVFull[, Scolnames]
colnames(T1)=c("CHROM", "START", "END", "CN", "Pheno", "ACMG_class","Location", "Pop.Freq", "Haploinsufficiency", "Triplosensitivity", "OMIM", "OMIM_morbid", "Reg.Element.Genes", "Nearest.Splice.Site.bp", "TAD", "BenignSource", "PathogenicSource", "gnomad.pLI", "exac.pLI", "GenesOfInterest", "Cosmic", "Pathways", "GTex", "Gene_name")

T1$CHROM=factor(T1$CHROM, levels=c(1:22, "X", "Y"))
T1=T1[order(T1$CHROM,T1$START), ]

############################################
# 4. Filter split length
############################################
Scolnames=c("SV_chrom", "SV_start", "SV_end","Gene_name", "Pathways","CN", "GenesOfInterest","Cosmic","Pheno", "ACMG_class", "Location", "AF", "HI", "TS","OMIM_phenotype", "OMIM_morb",
            "RE_gene", "Dist_nearest_SS", "TAD_coordinate", "Bsource", "Psource", "GnomAD_pLI", "ExAC_pLI", "GTex")       
# Filter based on whether there is a GeneofInterest or Cosmic related gene
idx=which(AnnotCNVSplit$GenesOfInterest!=""|AnnotCNVSplit$Cosmic!="")
T2=AnnotCNVSplit[idx, Scolnames]
colnames(T2)=c("CHROM","START","END", "Gene_name","Pathways","CN","PathwayRelated","Cosmic", "Pheno","ACMG_class","Location","Pop.Freq","Haploinsufficiency","Triplosensitivity",
               "OMIM","OMIM_morbid","Reg.Element.Genes","Nearest.Splice.Site.bp", "TAD","BenignSource","PathogenicSource","gnomad.pLI","exac.pLI","GTex")   
# also indicate whether to use a gene for a plot based on how long the CV is, or whether it crosses multiple intron/exon boundaries
SpanReg=strsplit(T2$Location,"-")
StartIdx=sapply(SpanReg, function(x) x[1])
EndIdx=sapply(SpanReg, function(x) x[2])
MultiInt=ifelse(StartIdx!=EndIdx, T, F)
MultiInt[grep("exon", StartIdx)]=T
MultiInt[grep("exon", EndIdx)]=T
T2$MultiSpan=MultiInt

T2$CHROM=factor(T2$CHROM, levels=c(1:22, "X", "Y"))
T2=T2[order(T2$CHROM,T2$START), ]

###############################################
# 5. Filter the ACMG genes
###############################################
Scolnames=c("SV_chrom", "SV_start", "SV_end", "CN", "Pheno", "ACMG_class", "AF", "Psource",
            "Gene_name","Cosmic", "Pathways","GTex","GenesOfInterest" )
filt2=which(AnnotCNVFull$ACMG_class>opts$ACMGcutoff)
T3=AnnotCNVFull[filt2, Scolnames]
colnames(T3)[c(1:3,7:8)]=c("CHROM", "START", "END", "Pop.Freq","PathogenicSource")

################################################
# 6. Write to file, or return the list
###############################################
  return(list(full=T1, split=T2, acmg=T3))
}