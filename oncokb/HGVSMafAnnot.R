#!/usr/bin/Rscript
"usage: \n HGVSMafAnnot.R [--maffile=<file> --outputfile=<string> --AAlist=<file> ]
\n options:\n --vcffile=<file> vcffile that has been annotated using vep.
\n --outputfile=<string> output directory
\n --sampleName=<string> If sample list has more than 1 sample, include a csv file with all identifiers, otherwise set as 'NULL' [default: NULL]
\n --AAlist=<file> File containing aminoacid 3 letter and 1 letter conversion " -> doc

library("docopt")
opts <- docopt(doc)

HGVSMafAnnot=function(maffile, outputfile, AAlist){
  AllData=read.delim(maffile, sep="\t", stringsAsFactors = F)
  AAlist=read.csv(AAlist)
  ## Create HGVSp short:
  print(colnames(AllData))
  print(match("HGVSp", colnames(AllData)))
  xa=strsplit(AllData[ ,match("HGVSp", colnames(AllData))], ":p\\.")
  ax2=sapply(xa, function(x) x[2])
  ax2=gsub("%3D", "=", ax2)
  ax3=ax2
  
  for (i in 1:nrow(AAlist)){
    ax3=gsub(AAlist$Three.letter.symbol[i], AAlist$One.letter.symbol[i], ax3)
  }
  
  print('Tidying HGVSp HGVSc')
  
  HGVSp_ENS=AllData[ ,match("HGVSp", colnames(AllData))]
  HGVSp=ax2
  ## change the sequence information
  HGVSp_Short=paste("p.", ax3 , sep="")
  HGVSp_Short[which(is.na(ax3))]=NA
  HGVSc_ENS=AllData[ ,match("HGVSc", colnames(AllData))]
  
  HGVSc=sapply(strsplit(AllData[, "HGVSc"], ":"), function(x) x[2])
  print('combine all information')
  Tx=AllData[ ,-match(c("HGVSc", "HGVSp"), colnames(AllData))]
  AllData=cbind(Tx,HGVSp_Short, HGVSp_ENS, HGVSp,HGVSc_ENS, HGVSc)
  colnames(AllData)[which(colnames(AllData)=="SYMBOL")]="Hugo_Symbol"
  print('write to file')
  write.table(AllData, file=outputfile, sep="\t", row.names=F, quote=F)
}

HGVSMafAnnot(opts$maffile, opts$outputfile, opts$AAlist)
