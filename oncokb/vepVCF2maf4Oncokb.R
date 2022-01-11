#!/usr/bin/Rscript
"usage: \n vepVCF2maf4Oncokb.R [--vcffile=<file> --outputfile=<string> --sampleName=<string> --AAlist=<file> --canonical=<boolean> ]
\n options:\n --vcffile=<file> vcffile that has been annotated using vep.
\n --outputfile=<string> output directory
\n --sampleName=<string> If sample list has more than 1 sample, include a csv file with all identifiers, otherwise set as 'NULL' [default: NULL]
\n --AAlist=<file> File containing aminoacid 3 letter and 1 letter conversion
\n --canonical=<boolean> Need to separate out canonical values if not run " -> doc

library("docopt")
opts <- docopt(doc)
#opts

## vcf to maf format
## following annotation with VEP104
## run with the options:
## --domains --af_gnomad --check_existing --fields CANONICAL,Consequence,Codons,Amino_acids,Gene,SYMBOL,Feature,EXON,PolyPhen,SIFT,Protein_position,BIOTYPE,Feature_type,cDNA_position,CDS_position,Existing_variation,DISTANCE,STRAND,CLIN_SIG,LoF_flags,LoF_filter,LoF,RadialSVM_score,RadialSVM_pred,LR_score,LR_pred,CADD_raw,CADD_phred,Reliability_index,DOMAINS,HGVSc,HGVSp,gnomAD_AF,PHENO

vcf2MAF=function(vcffile,outputfile,sampleName, AAlist, canonical=F){
  library(vcfR)
  library(matrixStats)
  AAlist=read.csv(AAlist)
  inputFile=read.vcfR(vcffile)
  vcfFix=getFIX(inputFile)
  # get depth
  dp1 <- extract.gt(inputFile, element='DP', as.numeric=TRUE)
  # if QUAL =0
    # get allelic depth: reference
  ad1 <- extract.gt(inputFile, element='AD', as.numeric=TRUE)
  af1=extract.gt(inputFile, element='AF', as.numeric=TRUE)
  adv <- dp1-ad1
  #vaf=adv/dp
  
  dp2<- extract.gt(inputFile, element='TIR', as.numeric=TRUE)
  ad2<- extract.gt(inputFile, element='TAR', as.numeric=TRUE)
  vaf2=ad2/(ad2+dp2)
  
  adv=ifelse(is.na(adv), ad2, adv)
  vaf=ifelse(is.na(af1), vaf2, af1)
  dp=ifelse(is.na(dp1), dp2, dp1)
  
  colnames(adv)=paste("nAlt", colnames(adv), sep="-")
  colnames(dp)=paste("nTot", colnames(dp), sep="-")
  colnames(vaf)=paste("VAF", colnames(vaf), sep="-")
  ## get the colnames for the CSQ upload
  a1=inputFile@meta
  x1=a1[grep("CSQ", a1)]
  SNames=unlist(strsplit( substr(x1, regexpr(  "Format: *",x1)+8, nchar(x1)-2), "\\|"))
  print('Extract information on CSQ')
  ann <- extract.info(inputFile, element='CSQ', as.numeric=F)
  annsplit=sapply(ann, function(x) strsplit(x, "\\||,"))
  print('Create new data table')
  ncounts=sapply(annsplit, length)
  
  print(SNames)
  
  if (length(unique(ncounts))>1){ 
    canonical = T }
  
  if (canonical==T){
    lxa=which(SNames=="CANONICAL")
    lxb=length(SNames)-lxa
    print(lxa)
    print(lxb)
    yesgrep=sapply(annsplit, function(x) grep("^YES$", x)) ## need to make sure protein sequences are not grepped
    tx=sapply(yesgrep, length)
    yesgrep2<-lapply(yesgrep, function(x) {ifelse(length(x)==0, x<-lxa, x); x})
    tx2=sapply(yesgrep2, length)
    csq3=sapply(1:length(annsplit), function(x) annsplit[[x]][(yesgrep2[[x]][1]-lxa+1): (yesgrep2[[x]][1]+lxb)]) ## issue line 932
    AllData=t(csq3)
    colnames(AllData)=SNames
  } else {
    names(annsplit)=paste("X",c(1:length(annsplit)))
    AllData=do.call(rbind, annsplit)
    colnames(AllData)=SNames[1:ncol(AllData)]    
  }

  ## Create HGVSp short:
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
   AllData=cbind(AllData[ ,-match(c("HGVSc", "HGVSp"), colnames(AllData))],HGVSp_Short, HGVSp_ENS, HGVSp,HGVSc_ENS, HGVSc)
   
  data.out=data.frame(Tumor_Sample_Barcode=sampleName, vcfFix, AllData, adv, dp, vaf)
  colnames(data.out)[which(colnames(data.out)=="SYMBOL")]="Hugo_Symbol"
  
  write.table(data.out, file=outputfile, sep="\t", row.names=F, quote=F)
}


#vcffile="~/Downloads/consensus_call-vep_MEL4_N.GRCh38_vep.vcf.gz" 
#vcf2MAF(vcffile, "~/Downloads/abra2_vc/vepcall/rerun_protein.maf", sampleName =  "ER099_MEL4N", AAlist="~/Documents/ER_pilot/annotations/AminoAcid_table.csv",protein=T,pfam="~/Documents/ER_pilot/annotations/Pfam-A.clans.33.1.tsv", pirsf = "~/Documents/ER_pilot/annotations/pirsfinfo.cat")

##print(opts)
vcf2MAF(opts$vcffile, opts$outputfile, opts$sampleName, opts$AAlist, opts$canonical)



