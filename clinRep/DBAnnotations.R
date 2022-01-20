#!/usr/bin/Rscript
"usage: \n DBAnnotations.R [--maffile=<file> --outputfile=<file> --cosmicMut=<file> --MSigDB=<file>]
\n options:\n --maffile=<file> maffile that has been annotated using vep.
\n --outputfile=<file> output file
\n --sampleName=<string> If sample list has more than 1 sample, include a csv file with all identifiers, otherwise set as 'NULL' [default: NULL]
\n --cosmicMut=<file> File containing list of Cosmic Mutations
\n --MSigDB=<file> File containing MsigDB gene set information " -> doc

library("docopt")
opts <- docopt(doc)

DBAnnotations=function(maffile, outputfile, cosmicMut, MSigDB){
  library(data.table)
  library(GSEABase)
  AllData=fread(maffile, sep="\t", stringsAsFactors = F, select=c("Hugo_Symbol",  "CHROM","POS", "HGVSp_Short","HGVSc"))
  print('annotate with MSigDB')

  PathInH=getGmt(con=MSigDB, geneIdType=SymbolIdentifier(),
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
  print('annotate with Cosmic')
  CosmicD=fread(cosmicMut, sep="\t", select=c("GENE_NAME","POS","Mutation AA","ONC_TSG","CGC_TIER","DISEASE", "CLINVAR_TRAIT","MUTATION_SIGNIFICANCE_TIER"))
  bbx1=paste(CosmicD$GENE_NAME, CosmicD$POS)
  bx1=paste(CosmicD$GENE_NAME, CosmicD$`Mutation AA`)
  mmx1=match(ax1, bx1)
  m1=match(aax1, bbx1)
  m1[which(!is.na(mmx1))]=mmx1[which(!is.na(mmx1))]
  CosmicD=CosmicD[m1,c("Mutation AA", "POS","ONC_TSG","CGC_TIER","DISEASE", "CLINVAR_TRAIT","MUTATION_SIGNIFICANCE_TIER") ]
  colnames(CosmicD)=paste("CMC", colnames(CosmicD), sep=".")
  AllData=cbind(AllData, CosmicD)
  print('write to file')
  write.table(AllData[ ,-c(1:4)], file=outputfile, sep="\t", quote=F,row.names = F)
}

DBAnnotations(opts$maffile,opts$outputfile, opts$cosmicMut, opts$MSigDB)
