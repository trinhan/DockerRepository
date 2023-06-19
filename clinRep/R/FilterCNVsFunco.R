#!/usr/bin/Rscript
"usage: \n FilterCNVs.R [--tsv=<file> --outputname=<string> --CNlow=<int> --CNhigh=<int> --ACMGcutoff=<int>]
\n options:
\n --tsv=<file> tsv from funcotator
\n --outputname=<string> output string
\n --CNlow=<int> [default: 1] is this in germline mode
\n --CNhigh=<int> [default: 3] is this in germline mode
\n --ACMGcutoff=<int> [default: 4] for germline mode " -> doc

library("docopt", quietly = T)
opts <- docopt(doc, help=TRUE, version='1.0.0')

#FilterCNVs=function(tsv, CNlow, CNhigh,ACMGcutoff ){
############################################
#1. Load in all required libraries silently
############################################

suppressMessages(library(data.table, quietly = T))
suppressMessages(library(matrixStats, quietly = T))

############################################
#2. Read in data file, set the CN thresholds and separate into split and full
############################################

CNV=read.delim(opts$tsv)
CNV2=CNV[which(CNV$Segment_Mean<opts$CNlow |CNV$Segment_Mean>opts$CNhigh), ]
sprintf('There are %s rows in data after filter less than %s and higher than %s ', nrow(CNV2), opts$CNlow, opts$CNhigh)
CNVFull=CNV2
# create a cnv split mode
CNVSplit=CNV2[CNV2$Annotation_mode=="split", ]


###########################################
# 3. Create a summary table
###########################################
message('creating summary table')
gain=which(CNVFull$SV_type=="DUP")
homloss=which(CNVFull$SV_type=="DEL")

VariantSummary=c("Number of DUP/GAINs ", "Number of DELs")

Nx=data.frame(Ngain=length(gain), Ndel=length(homloss))
SummTable=cbind(VariantSummary, t(Nx))

############################################
#2B. Prepare the dictionary for renaming data
############################################
#if (opts$ColsID=="NULL"){
  SIDCol=read.csv("./annotFiles/SV_column_IDs.csv")
#} else {
#  SIDCol=read.csv(opts$ColsID)
#}


############################################
#3. Filter full length
############################################
 scolidx=which(SIDCol$Table%in%c("ALL", "FULL"))
 m1=match(SIDCol$Variable[scolidx], colnames(CNVFull))
 T1=CNVFull[, na.omit(m1)]
 colnames(T1)=SIDCol$ReName[scolidx[which(!is.na(m1))]]
# ## Reannotate based on chromosome
 T1$CHROM=factor(T1$CHROM, levels=c(1:22, "X", "Y"))
 T1=T1[order(T1$CHROM,T1$START), ]

# Scolnames=c("SV_chrom", "SV_start", "SV_end", "CN", "Pheno", "ACMG_class", "Location", "AF", "HI", "TS","OMIM_phenotype", "OMIM_morb",
#             "RE_gene", "Dist_nearest_SS", "TAD_coordinate",
#             "Bsource", "Psource", "GnomAD_pLI", "ExAC_pLI","GenesOfInterest","Cosmic", "Pathways","GTex", "Gene_name")
# # Filter based on whether there is a GeneofInterest or Cosmic related gene
# # Add a tag to indicate the genes of interest
# T1=CNVFull[, Scolnames]
# colnames(T1)=c("CHROM", "START", "END", "CN", "Pheno", "ACMG_class","Location", "Pop.Freq", "Haploinsufficiency", "Triplosensitivity", "OMIM", "OMIM_morbid", "Reg.Element.Genes", "Nearest.Splice.Site.bp", "TAD", "BenignSource", "PathogenicSource", "gnomad.pLI", "exac.pLI", "GenesOfInterest", "Cosmic", "Pathways", "GTex", "Gene_name")
# 
# T1$CHROM=factor(T1$CHROM, levels=c(1:22, "X", "Y"))
# T1=T1[order(T1$CHROM,T1$START), ]

sprintf('Full length CNVs associated with pathway/Cosmic: %s entries', nrow(T1))


############################################
# 4. Filter split length
############################################
# Scolnames=c("SV_chrom", "SV_start", "SV_end","Gene_name", "Pathways","CN", "GenesOfInterest","Cosmic","Pheno", "ACMG_class", "Location", "AF", "HI", "TS","OMIM_phenotype", "OMIM_morb",
#             "RE_gene", "Dist_nearest_SS", "TAD_coordinate", "Bsource", "Psource", "GnomAD_pLI", "ExAC_pLI", "GTex")       
# # Filter based on whether there is a GeneofInterest or Cosmic related gene
# idx=which(CNVSplit$GenesOfInterest!=""|CNVSplit$Cosmic!="")
# T2=CNVSplit[idx, Scolnames]
# colnames(T2)=c("CHROM","START","END", "Gene_name","Pathways","CN","PathwayRelated","Cosmic", "Pheno","ACMG_class","Location","Pop.Freq","Haploinsufficiency","Triplosensitivity",
#                "OMIM","OMIM_morbid","Reg.Element.Genes","Nearest.Splice.Site.bp", "TAD","BenignSource","PathogenicSource","gnomad.pLI","exac.pLI","GTex")   
# # also indicate whether to use a gene for a plot based on how long the CV is, or whether it crosses multiple intron/exon boundaries
# SpanReg=strsplit(T2$Location,"-")
# StartIdx=sapply(SpanReg, function(x) x[1])
# EndIdx=sapply(SpanReg, function(x) x[2])
# MultiInt=ifelse(StartIdx!=EndIdx, T, F)
# MultiInt[grep("exon", StartIdx)]=T
# MultiInt[grep("exon", EndIdx)]=T
# T2$MultiSpan=MultiInt
# 
# T2$CHROM=factor(T2$CHROM, levels=c(1:22, "X", "Y"))
# T2=T2[order(T2$CHROM,T2$START), ]


scolidx=which(SIDCol$Table%in%c("ALL", "SPLIT"))
m1=match(SIDCol$Variable[scolidx], colnames(CNVSplit))
S2=CNVSplit[, na.omit(m1)]
colnames(S2)=SIDCol$ReName[scolidx[which(!is.na(m1))]]
# extract the specific genes of interest
idx=which(CNVSplit$GenesOfInterest!=""|CNVSplit$Cosmic!="")
try(if(length(idx)==0) warning("No genes of interest found"))
S2=S2[idx, ]

if(length(idx)>0){
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
}

sprintf('Specific genes associated with pathway/Cosmic: %s entries', nrow(S2))

###############################################
# 5. Filter the ACMG genes
###############################################
Scolnames=c("SV_chrom", "SV_start", "SV_end", "CN", "Pheno", "ACMG_class", "AF", "Psource",
            "Gene_name","Cosmic", "Pathways","GTex","GenesOfInterest" )
filt2=which(CNVFull$ACMG_class>opts$ACMGcutoff)
T3=CNVFull[filt2, Scolnames]
colnames(T3)[c(1:3,7:8)]=c("CHROM", "START", "END", "Pop.Freq","PathogenicSource")

sprintf('Specific genes associated with pathway/Cosmic: %s entries', nrow(T3))


################################################
# 6. Write to file, or return the list
###############################################
################################################
# 6. Write to file, or return the list
###############################################
message('write files to output')
##write.table(SummTable, file=paste(opts$outputname, ".CNV.SummaryTable.txt", sep=""), sep = "\t", row.names = F,  quote = F)
write.table(T1, file=paste(opts$outputname, ".CNV.full.filt.maf", sep=""), sep = "\t", row.names = F,  quote = F)
write.table(T2, file=paste(opts$outputname, ".CNV.split.filt.maf", sep=""), sep = "\t", row.names = F,  quote = F)
write.table(T3, file=paste(opts$outputname, ".CNV.acmg.filt.maf", sep=""), sep = "\t", row.names = F,  quote = F)

