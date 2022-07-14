#!/usr/bin/Rscript
"usage: \n SummarizeAnnotSV.R [--AnnotSVtsv=<file> --outputname=<string> --germline=<boolean> --MSigDB=<file> --GTex=<file> --CosmicList=<file> --AddList=<file> --pathwayList=<string> --ACMGCutoff=<integer> --Tissue=<string> --CNV=<boolean> --SRfilter=<integer> --PRfilter=<integer> --PASSfilt=<boolean>]
\n options:\n --AnnotSVtsv=<file> tsv from AnnotSV.
\n --outputname=<string> output string
\n --germline=<boolean> [default: TRUE] is this in germline mode
\n --MSigDB=<file> File containing MsigDB gene set information
\n --GTex=<file> GTex Expression Data
\n --CosmicList=<file> File of Cosmic Related Genes
\n --AddList=<file> List of additionalgene lists [default=NULL]
\n --pathwayList=<string>
\n --ACMGCutoff=<integer> Cut off value to assess these variants. Anything below this will not be further annotated [default:4]
\n --Tissue=<string> Tissue type 
\n --CNV=<boolean> Whether the input is a CNV or SV
\n --SRfilter=<integer> [default: 1] Minimum number of SR to support SV
\n --PRfilter=<integer> [default: 1] Minimum number of PR to support SV
\n --PASSfilt=<boolean> filter variants with PASS tag [default=TRUE] " -> doc

## Include the pathway List!

library("docopt", quietly = T)
opts <- docopt(doc)

###SummarizeAnnotSV=function(AnnotSVtsv, outputname, germline=T, MSigDB, GTex,CosmicList, GeneList=NULL, ACMGCutoff=4, Tissue="skin"){
  
  library(data.table, quietly = T)
  library(dplyr, quietly = T)
  library(GSEABase, quietly = T)
  library(matrixStats, quietly = T)
  InputData=read.delim(opts$AnnotSVtsv)
  
  if (opts$PASSfilt){
    InputData=InputData[which(InputData$FILTER=="PASS"| InputData$FILTER=="."), ]
  }
  
  print('Remove blacklisted regions from ENCODE')
  rmx=which(InputData$ENCODE_blacklist_characteristics_left=="" & InputData$ENCODE_blacklist_characteristics_right=="")
  InputData=InputData[rmx, ]
  print('Tidy pheno-geno data')
  Psource=ifelse(InputData$SV_type=="DEL", InputData$P_loss_source,ifelse(InputData$SV_type=="DUP", InputData$P_gain_source, InputData$P_ins_source))
  Bsource=ifelse(InputData$SV_type=="DEL", InputData$B_loss_source,
                 ifelse(InputData$SV_type=="DUP", InputData$B_gain_source,ifelse(InputData$SV_type=="INV", InputData$B_inv_source, InputData$B_ins_source)))
  BAF=ifelse(InputData$SV_type=="DEL",  InputData$B_loss_AFmax,
             ifelse(InputData$SV_type=="DUP", InputData$B_gain_AFmax, 
                    ifelse(InputData$SV_type=="INV", InputData$B_inv_AFmax, InputData$B_ins_AFmax)))
  
  ## Extract the genotype and copy number
  sampleName=opts$outputName ##data$sampleIn 
  temp=strsplit(as.character(InputData[ ,15 ]), ":")
  GT=sapply(temp, function(x) x[1])

  Pheno=ifelse(InputData$SV_type=="DEL", InputData$P_loss_phen, ifelse(InputData$SV_type=="DUP", InputData$P_gain_phen, 
                                                                       ifelse(InputData$SV_type=="INS", InputData$P_ins_phen, InputData$P_snvindel_phen)))
  
  InputData$GT=GT
  InputData$Pheno=Pheno
  InputData$AF=BAF
  InputData$Psource=Psource
  InputData$Bsource=Bsource
  Genes=strsplit(InputData$Gene_name, ";")
  ### Get the information on ACMG genes
  filtVals=which(InputData$ACMG_class>=opts$ACMGCutoff | InputData$TAD_coordinate!="")
  
  print('Annotate implicated genes with pathway data')  
  PathInH=getGmt(con=opts$MSigDB, geneIdType=SymbolIdentifier(),
                 collectionType=BroadCollection(category="h"))
  Mx2=geneIds(PathInH)
  GS=unique(unlist(Genes[filtVals]))
  GS2=paste("^", GS, "$", sep="")
  Ids2=names(PathInH)
  Ids2=gsub("HALLMARK_", "", Ids2)
  names(Mx2)=Ids2
  Mx3=stack(Mx2)
  Tx2=sapply(GS, function(x) as.character(Mx3$ind[which(Mx3$values==x)]))
  Nm2=sapply(Genes[filtVals], function(x) match(x, names(Tx2)))
  Nm3=sapply(Nm2, function(x) paste(unique(as.character(Tx2[x])), collapse=" "))
  Nm3=gsub( "character\\(0\\)", "",Nm3)
  Nm3=gsub("c\\(\"", "", Nm3)
  Nm3=gsub("\", \"", " ", Nm3)
  Nm3=gsub("\")", "", Nm3)
  
  InputData$Pathways=NA
  InputData$Pathways[filtVals]=Nm3
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
  ## also include information from the GeneList
  if (!is.null(opts$AddList)){
    print('Annotate with input gene list')
    GL=read.delim(opts$AddList, header=F)
    GL=GL[ ,1]
    Nx2=sapply(Genes, function(x) paste(x[which(x%in%GL)], collapse=", "))
    InputData$GenesOfInterest=Nx2
  }else{
    InputData$GenesOfInterest=NA
  }
  
  sprintf('Find genes involved in %s pathway', opts$pathwayList)
  ## Pathways of Interest
  #PWtable=read.csv("~/Documents/ER_pilot/annotations//PathwayList.csv")
  PWtable=read.csv("/annotFiles/PathwayList.csv")
  nx=which(colnames(PWtable)==opts$pathwayList)
  nx=setdiff(as.vector(PWtable[ ,nx]), "")
  ex1=unique(unlist(sapply(nx, function(x) grep(x, InputData$Pathways))))
  InputData$GenesOfInterest[ex1]=ifelse(is.na(InputData$GenesOfInterest[ex1]), opts$pathwayList, paste(InputData$GenesOfInterest[ex1], opts$pathwayList))
  
  print('Add Cosmic Data')
  ## CosmicData
  CData=read.csv(opts$CosmicList)
  Nx2=sapply(Genes, function(x) paste(x[which(x%in%CData$Gene.Symbol)], collapse=", "))
  InputData$Cosmic=Nx2
  
  print('Additional SV or CNV specific filter')
  ## Additional information based on AnnotSV output
  if (opts$CNV){
    CN=sapply(temp, function(x) x[2])
    InputData$CN=CN
  } else {
    PR=sapply(temp, function(x) x[5])
    PR2=strsplit(PR, ",")
    SR=sapply(temp, function(x) x[6])
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
    tx1=which(PRcounts>opts$PRfilter & SRcounts>opts$SRfilter)
    InputData=InputData[tx1, ]
  }
  
  
  print('Finished! Write out formated output table')
  if (opts$CNV){
    write.table(InputData, file=paste(opts$outputname, ".CNV.formated.tsv", sep=""), sep="\t", row.names = F, quote = F)
  } else {
    InputDataFull=InputData[InputData$Annotation_mode=="full", ]
    InputDataSplit=InputData[InputData$Annotation_mode=="split", ]
    write.table(InputDataFull, file=paste(opts$outputname, ".SV.formated.tsv", sep=""), sep="\t", row.names = F, quote = F)
    write.table(InputDataSplit, file=paste(opts$outputname, ".SV.split.formated.tsv", sep=""), sep="\t", row.names = F, quote = F)
  }
#}

##SummarizeAnnotSV(opts$AnnotSVtsv, opts$outputname, opts$germline, opts$MSigDB, opts$GTex,opts$CosmicList, opts$GeneList, opts$ACMGCutoff,
##                 opts$Tissue)