#!/usr/bin/Rscript
"
This script converts a vep.vcf output file into a .maf file

Usage: \n 
vepVCF2maf.R --vcffile=<file> --outputfile=<string> --sampleName=<string> [--canonical=<boolean> --runMode=<string> --AAlist=<file> --minCallers=<string>]
\n options:\n --vcffile=<file> vcffile that has been annotated using vep.
\n --outputfile=<string> output directory
\n --sampleName=<string> If sample list has more than 1 sample, include a csv file with all identifiers, otherwise set as 'NULL' [default: NULL]
\n --canonical=<boolean> Need to identify the canonical transcript if not run [default: T]
\n --runMode=<string> Tumour or germline mode? [default: 'Tumour']
\n --AAlist=<file> File containing aminoacid 3 letter and 1 letter conversion 
" -> doc

library("docopt")
opts <- docopt(doc)
print(opts)


############################################
#1. Load in all required libraries silently
############################################
suppressMessages(library(vcfR))
suppressMessages(library(matrixStats))

#####################################################################
#2. Read in the vcf file
#####################################################################  

print('Running vep VCF to MAF')
inputFile=read.vcfR(opts$vcffile)
vcfFix=getFIX(inputFile)

#####################################################################
#3. Get depth and read count information
#####################################################################  

if (opts$runMode=="Tumour"){
  print('Get depth info for tumour')
  # SNV-extract total depth
  dp1 <- extract.gt(inputFile, element='DP', as.numeric=TRUE)
  # SNV-Extract the alt counts
  ad2 <- extract.gt(inputFile,  element='AD',as.numeric=T)
  adv1 <- dp1-ad2
  # SNV-calclate af
  af1=extract.gt(inputFile, element='AF', as.numeric=TRUE)
  # extract information on Strelka runs here
  
  ad2<- extract.gt(inputFile, element='TIR', as.numeric=TRUE)
  adv2<- extract.gt(inputFile, element='TAR', as.numeric=TRUE)
  vaf2=adv2/(adv2+ad2)
  print('calculate for strelka runs')
  dp2=adv2+ad2
  
  # calculate nAlt reads
  adv=ifelse(is.na(adv1), adv2, adv1)
  # calculate vaf: take the default AF tag or compute it for strelka
  vaf=ifelse(is.na(af1), vaf2, af1)
  # calculate depth - take either default or that obtained from strelka
  dp=ifelse(is.na(dp1), dp2, dp1)
  # obtain the genotype
  gt<-extract.gt(inputFile, element='GT')
} else if (opts$runMode=="Germline") {
  print('Get depth info for normal')
  # extract the AD tag
  adv <- extract.gt(inputFile, element='AD', as.numeric=F)
  t2=colnames(adv)
  # split into ref vs alt counts
  adA = strsplit(adv, ",")
  fwdA=sapply(adA, function(x) x[2])
  revA=sapply(adA, function(x) x[1])
  adv=matrix(as.numeric(fwdA), ncol=ncol(adv))
  colnames(adv)=t2
  adv2=matrix(as.numeric(revA), ncol=ncol(adv))
  colnames(adv2)=t2
  # compute depth and vaf (some callers do not have this info)
  dp=adv+adv2
  vaf=adv/dp
  gt<-extract.gt(inputFile, element='GT')
}

colnames(adv)=paste("nAlt", colnames(adv), sep="-")
colnames(dp)=paste("nTot", colnames(dp), sep="-")
colnames(vaf)=paste("VAF", colnames(vaf), sep="-")
colnames(gt)=paste("GT", colnames(gt), sep="-")

VAFd=cbind(adv, dp, vaf, gt)

#####################################################################
#4. Preparation to extract VEP annotations
#####################################################################  

LxB=nrow(vcfFix)
sprintf('Get CSQ colnames... There are %s variants', LxB)

# get the names of the vep columns
a1=inputFile@meta
x1=a1[grep("CSQ", a1)]
SNames=unlist(strsplit( substr(x1, regexpr(  "Format: *",x1)+8, nchar(x1)-2), "\\|"))
# read in AA amino acid file
AAlist=read.csv(opts$AAlist)

# Note that reading the whole vcf file and storing in memory can take a long time.
# Here, we try to parallelise this by reading in 50000 rows at a time

# function to define how many for loops are required for reading.
seqlast <- function (from, to, by) 
{
  vec <- do.call(what = seq, args = list(from, to, by))
  if ( tail(vec, 1) != to ) {
    return(c(vec, to))
  } else {
    return(vec)
  }
}

Ntries <- seqlast(from=0, to=LxB, by = 50000)

#####################################################################
# 5.Extract CSQ data and save to file 50000 lines at a time 
#####################################################################  
print('Create new data table')

# start the loop
for (i in 1:(length(Ntries)-1)){
  #Extract the CSQ information and save as a list
  ann <- extract.info(inputFile[(Ntries[i]+1):Ntries[i+1]], element='CSQ', as.numeric=F)
  annsplit<-sapply(ann, function(x) strsplit(x, "\\||,"))
  ncounts<- sapply(annsplit, length)
  if (opts$canonical){
    # run this section if we only want the canonical sample
    # Note this section only remains due to previous coding decisions. Can decide to delete
    print('find the canonical gene.')
    # figure out which column header is the canonical sample 
    lxa=which(SNames=="CANONICAL")
    lxb=length(SNames)-lxa
    # grep "YES" to find the index of the canonical sample
    yesgrep=sapply(annsplit, function(x) grep("^YES$", x))
    # for the entries within a canonical sample, assign the canonical sample the default column number (first transcript entry)
    yesgrep2<-lapply(yesgrep, function(x) {ifelse(length(x)==0, x<-lxa, x); x})
    csq3=sapply(1:length(annsplit), function(x) annsplit[[x]][(yesgrep2[[x]][1]-lxa+1): (yesgrep2[[x]][1]+lxb)])
    AllData=t(csq3)
    colnames(AllData)=SNames
  } else if (length(unique(ncounts))>1) {
    # run this section if do not care about the canonical sample, but the number of CSQ entries is high. 
    # Note this section only remains due to previous coding decisions. Can decide to delete
    print('Do not use selection criteria')
    csq3=sapply(1:length(annsplit), function(x) annsplit[[x]][1:length(SNames)]) ## issue line 932
    AllData=t(csq3)
    colnames(AllData)=SNames
  }  else {
    print ('No canonical search required')
    # do a default renaming of columns, otherwise it will take the default which is very long
    names(annsplit)=paste("X",c(1:length(annsplit)))
    AllData=do.call(rbind, annsplit)
    # rename the columns
    colnames(AllData)=SNames[1:ncol(AllData)]    
  }
  print('Merge CSQ data')
  AllData=cbind(Tumor_Sample_Barcode=opts$sampleName, vcfFix[(Ntries[i]+1):Ntries[i+1], ], AllData, VAFd[(Ntries[i]+1):Ntries[i+1], ])
  
  print('Tidying HGVSp HGVSc')
  ## Create HGVSp short:
  xa=strsplit(AllData[ ,match("HGVSp", colnames(AllData))], ":p\\.")
  ax2=sapply(xa, function(x) x[2])
  ax2=gsub("%3D", "=", ax2)
  ax3=ax2 # keep the original version for gsubing later
  
  for (i in 1:nrow(AAlist)){
    ax3=gsub(AAlist$Three.letter.symbol[i], AAlist$One.letter.symbol[i], ax3)
  }
  
  # define additional columns:
  # a. HGVSp ensembl
  HGVSp_ENS=AllData[ ,match("HGVSp", colnames(AllData))]
  # b. HGVSp original
  HGVSp=ax2
  # c. HGVSp short - needed for oncokb look-up
  ## change the sequence information
  HGVSp_Short=paste("p.", ax3 , sep="")
  HGVSp_Short[which(is.na(ax3))]=NA
  # d. HGVSc ensembl
  HGVSc_ENS=AllData[ ,match("HGVSc", colnames(AllData))]
  HGVSc=sapply(strsplit(AllData[, "HGVSc"], ":"), function(x) x[2])
  print('combine all information and white to file')
  Tx=AllData[ ,-match(c("HGVSc", "HGVSp"), colnames(AllData))]
  AllData=cbind(Tx,HGVSp_Short, HGVSp_ENS, HGVSp,HGVSc_ENS, HGVSc)
  colnames(AllData)[which(colnames(AllData)=="SYMBOL")]="Hugo_Symbol"
  sprintf('write to file. Table dimensions is %s rows by %s columns', nrow(AllData), ncol(AllData))
  write.table(AllData, file=opts$outputfile, sep="\t", row.names=F, quote=F, append = T,
              col.names=!file.exists(opts$outputfile))
}
