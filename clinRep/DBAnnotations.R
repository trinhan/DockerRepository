#!/usr/bin/Rscript
"usage: \n DBAnnotations.R [--maffile=<file> --outputfile=<file> --cosmicMut=<file> --cosmicGenes=<file> --MSigDB=<file> --pfam=<file> --pirsf=<file>]
\n options:\n --maffile=<file> maffile that has been annotated using vep.
\n --outputfile=<file> output file
\n --sampleName=<string> If sample list has more than 1 sample, include a csv file with all identifiers, otherwise set as 'NULL' [default: NULL]
\n --cosmicMut=<file> File containing list of Cosmic Mutations
\n --cosmicGenes=<file> File of Cosmic Cancer Genes
\n --MSigDB=<file> File containing MsigDB gene set information 
\n --pfam=<file> File containing pfam annotations [default: NULL]
\n --pirsf=<file> File containing psirf annotations [default: NULL] " -> doc

library("docopt", quietly = T)
opts <- docopt(doc)

  library(data.table, quietly = T)
  library(GSEABase, quietly = T)
 library(dplyr, quietly = T)
  
  print('reading in maffile')
  
  AllData=fread(opts$maffile, sep="\t", stringsAsFactors = F, select=c("Hugo_Symbol",  
                          "CHROM","POS", "REF", "ALT", "HGVSp_Short","HGVSc", "DOMAINS"))
  print('annotate with MSigDB')

  PathInH=getGmt(con=opts$MSigDB, geneIdType=SymbolIdentifier(),
                 collectionType=BroadCollection(category="h"))
  Mx2=geneIds(PathInH)
  GS=unique(AllData$Hugo_Symbol)
  GS2=paste("^", GS, "$", sep="")
  Ids2=names(PathInH)
  Ids2=gsub("HALLMARK_", "", Ids2)
  names(Mx2)=Ids2
  Mx3=stack(Mx2)
  Tx2=sapply(GS, function(x) paste(Mx3$ind[which(Mx3$values==x)], collapse=" "))
  Nm2=match(AllData$Hugo_Symbol, GS)
  AllData$HallmarkPathways=Tx2[Nm2]
  
  aax1=paste(AllData$Hugo_Symbol, AllData$POS)
  ax1=paste(AllData$Hugo_Symbol, AllData$HGVSp_Short)
  rm(PathInH)
  print('annotate with Cosmic')
  CosmicD=fread(opts$cosmicMut, sep="\t", select=c("GENE_NAME","POS","Mutation AA","ONC_TSG","CGC_TIER","DISEASE", "CLINVAR_TRAIT","MUTATION_SIGNIFICANCE_TIER"))
  bbx1=paste(CosmicD$GENE_NAME, CosmicD$POS)
  bx1=paste(CosmicD$GENE_NAME, CosmicD$`Mutation AA`)
  mmx1=match(ax1, bx1)
  m1=match(aax1, bbx1)
  m1[which(!is.na(mmx1))]=mmx1[which(!is.na(mmx1))]
  
  CosmicD=CosmicD[m1,c("Mutation AA", "POS","ONC_TSG","CGC_TIER","DISEASE", "CLINVAR_TRAIT","MUTATION_SIGNIFICANCE_TIER") ]
  colnames(CosmicD)=paste("CMC", colnames(CosmicD), sep=".")
  AllData=cbind(AllData, CosmicD)
  # another comment: find genes which are oncogenes, but not exact mut site
  rm(CosmicD)
  
  print('check Cosmic Genes')
  Cosmic2=read.csv(opts$cosmicGenes)
  tx2=match(AllData$Hugo_Symbol, Cosmic2$Gene.Symbol)
  #head(Cosmic2)
  #length(na.omit(tx2))
  AllData$CMC.Cancer_Gene_Tier=NA
  AllData$CMC.Cancer_Gene_Tier[which(!is.na(tx2))]=Cosmic2$Tier[na.omit(tx2)]
  tx3=match(AllData$Hugo_Symbol, Cosmic2$Gene.Symbol[which(Cosmic2$Hallmark=="Yes")])
  AllData$CMC.Cancer_Gene_Tier[which(!is.na(tx3))]=paste("Hallmark", AllData$CMC.Cancer_Gene_Tier[which(!is.na(tx3))])
  
  print('Add protein domain annotations from pfam and psird')
  
  pfamI=read.delim(opts$pfam, header=F, sep="\t")

  pirsfinfo=read.delim(opts$pirsf, sep=")", header=F, quote="")
  #head(pirsfinfo)
  #class(pirsfinfo)
  pirsfinfo$V3=NA
  pirsfinfo$V3=substr(pirsfinfo$V1, 2, 12)
  
  ##pirsfinfo$V1=gsub(">", "", pirsfinfo$V1)
  print('add pfam')
  n2=regmatches( AllData$DOMAINS, 
                 regexec("PF[0-9]{5}", AllData$DOMAINS))
  n2pfam=unlist({n2[sapply(n2, length)==0] <- NA; n2})
  pfamO=as.character(pfamI$V5[match(n2pfam, pfamI$V1)])
  
  print('add pirsf')
  n2=regmatches( AllData$DOMAINS, 
                 regexec("PIRSF[0-9]{6}", AllData$DOMAINS))
  n2=unlist({n2[sapply(n2, length)==0] <- NA; n2})
  pirsfO=as.character(pirsfinfo$V2[match(n2,pirsfinfo$V3)])
  
  AllData=cbind(AllData, pfamid=n2pfam, pfam=pfamO, pirsf=pirsfO, pirsfid=n2)
 ## write.table(AllData, file=opts$outputfile, sep="\t", row.names=F, quote=F)
  print('write to file')
  write.table(AllData[ ,-c(1:4)], file=opts$outputfile, sep="\t", quote=F,row.names = F)
#}

