#!/usr/bin/Rscript
"usage: \n SummarizeAnnotSV.R [--AnnotSVtsv=<file> --outputname=<string> --germline=<boolean> --MSigDB=<file> --GTex=<file> --CosmicList=<file> --AddList=<file> --pathwayTerm=<string> --pathwayList=<file> --ACMGCutoff=<integer> --Tissue=<string> --CNV=<boolean> --SRfilter=<integer> --PRfilter=<integer> --PASSfilt=<boolean>]
\n options:\n --AnnotSVtsv=<file> tsv from AnnotSV.
\n --outputname=<string> output string
\n --germline=<boolean> [default: TRUE] is this in germline mode
\n --MSigDB=<file> File containing MsigDB gene set information
\n --GTex=<file> GTex Expression Data
\n --CosmicList=<file> File of Cosmic Related Genes
\n --AddList=<file> List of additionalgene lists [default=NULL]
\n --pathwayTerm=<string>
\n --pathwayList=<file>
\n --ACMGCutoff=<integer> Cut off value to assess these variants. Anything below this will not be further annotated [default:3]
\n --Tissue=<string> Tissue type 
\n --CNV=<boolean> Whether the input is a CNV or SV
\n --SRfilter=<integer> [default: 1] Minimum number of SR to support SV
\n --PRfilter=<integer> [default: 1] Minimum number of PR to support SV
\n --PASSfilt=<boolean> filter variants with PASS tag [default=TRUE] " -> doc

library("docopt", quietly = T)
opts <- docopt(doc)

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

InputData=read.delim(opts$AnnotSVtsv)
print(paste0("Number of rows:", nrow(InputData)))
  
if (opts$PASSfilt){
  InputData=InputData[which(InputData$FILTER=="PASS"| InputData$FILTER=="."), ]
  print(paste0("Number of rows PASS filter:", nrow(InputData)))
}

nrowCheck=(nrow(InputData)<=1)

# Put in a check to stop the script if everything is filtered out
{
  if (nrowCheck) {stop("There are insufficient lines in this file")}
}

############################################
#3. Filter out Encode Blacklisted regions
############################################

# Remove all regions which are blacklisted and have annotations
print('Remove blacklisted regions from ENCODE')
keepidx=which((InputData$ENCODE_blacklist_characteristics_left==""|is.na(InputData$ENCODE_blacklist_characteristics_left)) & 
              (InputData$ENCODE_blacklist_characteristics_right==""|is.na(InputData$ENCODE_blacklist_characteristics_right)))
InputData=InputData[keepidx, ]
print(paste0("Number of rows after Encode review:", nrow(InputData)))

############################################
#4. Filter out regions which do not fulfil the ACMG score cut-off.
#   If there is no ACMG cut-off, make sure it overlaps a TAD region
############################################

filtVals=which(InputData$ACMG_class>=opts$ACMGCutoff | InputData$TAD_coordinate!="")

InputData=InputData[filtVals, ]
print(paste0("Number of rows after ACMG review:", nrow(InputData)))

############################################
#5. Tidy up annotations on pathogenicity/benign
############################################

# Instead of having two columns with annotations on deletions and insertions, combine this information into 1 column based on the SV type
print('Tidy pheno-geno data')

## Tidy up the annotation source
# Psource - assumes there's just Deletions, Duplications and Insertions
Psource=ifelse(InputData$SV_type=="DEL", InputData$P_loss_source,ifelse(InputData$SV_type=="DUP", InputData$P_gain_source, InputData$P_ins_source))
# Bsource - assumes there's just Deletions, Duplications/Inserions and Inversions
Bsource=ifelse(InputData$SV_type=="DEL", InputData$B_loss_source,
                 ifelse(InputData$SV_type=="DUP", InputData$B_gain_source,ifelse(InputData$SV_type=="INV", InputData$B_inv_source, InputData$B_ins_source)))
# AF only avilable for Benign SVs
BAF=ifelse(InputData$SV_type=="DEL", InputData$B_loss_AFmax,ifelse(InputData$SV_type=="DUP", InputData$B_gain_AFmax, 
                    ifelse(InputData$SV_type=="INV", InputData$B_inv_AFmax, InputData$B_ins_AFmax)))
# Tidy-up pathogenic prediction
Pheno=ifelse(InputData$SV_type=="DEL", InputData$P_loss_phen, ifelse(InputData$SV_type=="DUP", InputData$P_gain_phen, 
                                                                     ifelse(InputData$SV_type=="INS", InputData$P_ins_phen, InputData$P_snvindel_phen)))
InputData$Pheno=Pheno
InputData$AF=BAF
InputData$Psource=Psource
InputData$Bsource=Bsource

############################################
#6. Extract the exact genes to search pathways - MSIGSB
############################################
print('Annotate implicated genes with pathway data')  

# Based on MSigDB
Genes=strsplit(InputData$Gene_name, ";")
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
if (!is.null(opts$AddList)){
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
#10. Extract genotypes for germline cases and CNs 
############################################

# Figure out how many samples are present
Nsamp=unlist(strsplit(InputData[ 1,"Samples_ID"], ","))
Nsamp2=Nsamp[length(Nsamp)]
Nsamp2=gsub("-", ".", Nsamp2)
print(Nsamp2)
  
temp=strsplit(as.character(InputData[ ,match(Nsamp2, colnames(InputData)) ]), ":")
head(temp)
  
if (opts$germline){
    GT=sapply(temp, function(x) x[1]) ## Check whether this is applicable for CNVs too?
    InputData$GT=GT
}

if (opts$CNV){
  print('CNV mode')
  CN=sapply(temp, function(x) x[2])
  InputData$CN=CN
}

############################################
#11. Extract PR and SR filters
############################################
print('Additional SV filter for PR and SR ')
  ## Additional information based on AnnotSV output
if (!opts$CNV){
    if (opts$germline){
      print('germline SV mode')
      PR=sapply(temp, function(x) x[5])
      SR=sapply(temp, function(x) x[6])
    }else{
      print('tumour SV mode')
      PR=sapply(temp, function(x) x[1])
      SR=sapply(temp, function(x) x[2])
    }  
    PR2=strsplit(PR, ",")
    SR2=strsplit(SR, ",")
    PRVAF=round(sapply(PR2, function(x) as.numeric(x[2])/(sum(as.numeric(x)))), 2)
    SRVAF=round(sapply(SR2, function(x) as.numeric(x[2])/(sum(as.numeric(x)))), 2)
    MeanVAF=round(sapply(1:length(PR2), function(x) sum(c(as.numeric(PR2[[x]][2]),as.numeric(SR2[[x]][2])), na.rm=T)/
                           sum(c(as.numeric(PR2[[x]]), as.numeric(SR2[[x]])), na.rm=T)), 2)
    InputData$PRVAF=PRVAF
    InputData$SRVAF=SRVAF
    InputData$MeanVAF=MeanVAF
    InputData$Depth=sapply(1:length(PRVAF), function(x) sum(c(as.numeric(PR2[[x]]), as.numeric(SR2[[x]])), na.rm = T))
    PRcounts=sapply(PR2, function(x) as.numeric(x[2]))
    SRcounts=sapply(SR2, function(x) as.numeric(x[2]))
    tx1=which(PRcounts>= as.numeric(opts$PRfilter) & SRcounts>=as.numeric(opts$SRfilter))
    print(max(na.omit(PRcounts)))
    print(max(na.omit(SRcounts)))
    InputData=InputData[tx1, ]
    print(paste0("Number of rows passing after SR and PR filters", nrow(InputData)))
}

  print('Finished! Write out formated output table')
if (opts$CNV){
    write.table(InputData, file=paste(opts$outputname, ".CNV.formated.tsv", sep=""), sep="\t", row.names = F, quote = F)
  } else {
    InputDataFull=InputData[InputData$Annotation_mode=="full", ]
    print(paste0("Number of SV rows in full:", nrow(InputDataFull)))
    InputDataSplit=InputData[InputData$Annotation_mode=="split", ]
    print(paste0("Number of SV rows in split:", nrow(InputDataSplit)))
    write.table(InputDataFull, file=paste(opts$outputname, ".SV.formated.tsv", sep=""), sep="\t", row.names = F, quote = F)
    write.table(InputDataSplit, file=paste(opts$outputname, ".SV.split.formated.tsv", sep=""), sep="\t", row.names = F, quote = F)
}
