#!/usr/bin/Rscript
"
This script converts a vep.vcf output file into a .maf file

Usage:
  vepVCF2maf.R --vcffile=<file> --outputfile=<string> --sampleName=<string> [--canonical=<boolean>] [--runMode=<string>] [--AAlist=<file>] [--minCallers=<string>]

Options:
  --vcffile=<file>       Vcffile that has been annotated using vep.
  --outputfile=<string>  Output directory.
  --sampleName=<string>  If sample list has more than 1 sample, include a csv file with all identifiers, otherwise set as 'NULL' [default: NULL].
  --canonical=<boolean>  Need to identify the canonical transcript if not run [default: T].
  --runMode=<string>     Tumour or germline mode? [default: 'Tumour'].
  --AAlist=<file>        File containing aminoacid 3 letter and 1 letter conversion.
  --minCallers=<string>  Minimum callers? [default: NULL].
" -> doc

suppressMessages(library("docopt", quietly = T))
suppressMessages(library("vcfR", quietly = T))
suppressMessages(library("matrixStats", quietly = T))

opts <- docopt(doc, help = TRUE, version = '1.0.0')
print(opts)

getDepthTumour<-function(inputFile){
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
  #calculate for strelka runs
  dp2=adv2+ad2
  
  # calculate nAlt reads
  adv=ifelse(is.na(adv1), adv2, adv1)
  # calculate vaf: take the default AF tag or compute it for strelka
  vaf=ifelse(is.na(af1), vaf2, af1)
  # calculate depth - take either default or that obtained from strelka
  dp=ifelse(is.na(dp1), dp2, dp1)
  # obtain the genotype
  gt<-extract.gt(inputFile, element='GT')
  colnames(adv)=paste("nAlt", colnames(adv), sep="-")
  colnames(dp)=paste("nTot", colnames(dp), sep="-")
  colnames(vaf)=paste("VAF", colnames(vaf), sep="-")
  colnames(gt)=paste("GT", colnames(gt), sep="-")
  
  VAFd<-cbind(adv, dp, vaf, gt)
}

getDepthGermline<-function(inputFile){
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
  colnames(adv)=paste("nAlt", colnames(adv), sep="-")
  colnames(dp)=paste("nTot", colnames(dp), sep="-")
  colnames(vaf)=paste("VAF", colnames(vaf), sep="-")
  colnames(gt)=paste("GT", colnames(gt), sep="-")
  
  VAFd=cbind(adv, dp, vaf, gt)
}

# Define function to determine number of loops for reading
seqlast <- function(from, to, by) {
  vec <- do.call(what = seq, args = list(from, to, by))
  if (tail(vec, 1) != to) {
    return(c(vec, to))
  } else {
    return(vec)
  }
}

# main function for extraction the CSQ columns
ExtractCSQ<-function(inputFile, SNames, canonical){
    ann <- extract.info(inputFile, element = 'CSQ', as.numeric = FALSE)
    annsplit <- sapply(ann, function(x) strsplit(x, "\\||,"))
    ncounts <- sapply(annsplit, length)
    
    if (canonical) {
      lxa <- which(SNames == "CANONICAL")
      lxb <- length(SNames) - lxa
      yesgrep <- sapply(annsplit, function(x) grep("^YES$", x))
      yesgrep2 <- lapply(yesgrep, function(x) {
        ifelse(length(x) == 0, x <- lxa, x)
        x
      })
      csq3 <- sapply(1:length(annsplit), function(x) annsplit[[x]][(yesgrep2[[x]][1] - lxa + 1):(yesgrep2[[x]][1] + lxb)])
      AllData <- t(csq3)
      colnames(AllData) <- SNames
    } else if (length(unique(ncounts)) > 1) {
      csq3 <- sapply(1:length(annsplit), function(x) annsplit[[x]][1:length(SNames)])
      AllData <- t(csq3)
      colnames(AllData) <- SNames
    } else {
      names(annsplit) <- paste("X", c(1:length(annsplit)))
      AllAllData <- do.call(rbind, annsplit)
      colnames(AllData) <- SNames[1:ncol(AllData)]
    }
      AllData
}

# Function to modify the Amino acid names

ChangeHGVSp<-function(AllData, AAlistData){
  xa <- strsplit(AllData[, match("HGVSp", colnames(AllData))], ":p\\.")
  ax2 <- sapply(xa, function(x) x[2])
  ax2 <- gsub("%3D", "=", ax2)
  ax3 <- ax2
  
  for (i in 1:nrow(AAlistData)) {
    ax3 <- gsub(AAlistData$Three.letter.symbol[i], AAlistData$One.letter.symbol[i], ax3)
  }
  
  HGVSp_ENS <- AllData[, match("HGVSp", colnames(AllData))]
  HGVSp <- ax2
  HGVSp_Short <- paste("p.", ax3, sep = "")
  HGVSp_Short[which(is.na(ax3))] <- NA
  HGVSc_ENS <- AllData[, match("HGVSc", colnames(AllData))]
  HGVSc <- sapply(strsplit(AllData[, "HGVSc"], ":"), function(x) x[2])
  
  Tx <- AllData[, -match(c("HGVSc", "HGVSp"), colnames(AllData))]
  AllData <- cbind(Tx, HGVSp_Short, HGVSp_ENS, HGVSp, HGVSc_ENS, HGVSc)
  colnames(AllData)[which(colnames(AllData) == "SYMBOL")] <- "Hugo_Symbol"
  AllData
}

# Function to convert vep.vcf to .maf file
convertVCFtoMAF <- function(vcffile, outputfile, sampleName, canonical = TRUE, runMode = "Tumour", AAlist = NULL, minCallers = NULL) {
  # Read the vcf file
  inputFile <- read.vcfR(vcffile)
  vcfFix <- getFIX(inputFile)
  
  # Get depth and read count information
  if (runMode == "Tumour") {
    VAFd <- getDepthTumour(inputFile)
  } else if (runMode == "Germline") {
    VAFd <- getDepthGermline(inputFile)
  }
  
  # Prepare to extract VEP annotations
  AAlistData <- read.csv(AAlist)
  
  # calculate the number of iterations needed
  LxB <- nrow(vcfFix)
  print(LxB)
  x1 <- inputFile@meta[grep("CSQ", inputFile@meta)]
  SNames <- unlist(strsplit(substr(x1, regexpr("Format: *", x1) + 8, nchar(x1) - 2), "\\|"))
  Ntries <- seqlast(from = 0, to = LxB, by = 50000)
  print(Ntries)
  
  # Create new data table
  AllData <- NULL
  # Start the reading loop
  for (i in 1:(length(Ntries) - 1)) {
    print(i)
    AllData<-ExtractCSQ(inputFile[ (Ntries[i]+1): Ntries[i+1],  ], SNames, canonical)
    head(AllData)
    AllData<-cbind(Tumor_Sample_Barcode = sampleName, vcfFix[(Ntries[i] + 1):Ntries[i + 1], ], AllData, VAFd[(Ntries[i] + 1):Ntries[i + 1], ])
    AllData<-ChangeHGVSp(AllData, AAlistData)
    write.table(AllData, file = outputfile, sep = "\t", row.names = FALSE, quote = FALSE, append = TRUE, col.names = !file.exists(outputfile))
    }
  print('VEPVCF2MAF COMPLETE') 
}

# Main function
convertVCFtoMAF(vcffile=opts$vcffile, outputfile=opts$outputfile, sampleName=opts$sampleName, 
                canonical=opts$canonical, runMode=opts$runMode, opts$AAlist, opts$minCallers)








      