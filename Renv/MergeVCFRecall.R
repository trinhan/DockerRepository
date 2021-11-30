#!/usr/bin/Rscript
# R script for altering the fields of a VCF file from judgement callers: blat, panels of normals (pon, gnomad) and abra2 indel realignment ()
"usage: \n MergeVCFRecall.R [--vcffile=<file> --outputstr=<string>  --genome=<string> --runMode=<string> --blat=<file> --pon=<file> --abra2=<file> --gnomad=<file> ]
\n options:\n --vcffile=<file> vcffile
\n --outputstr=<string> outputName
\n --genome=<string> reference genome: hg19 or hg38 [default: 'hg19']
\n --runMode=<string> run mode: TumOnly or Paired [default: 'Paired']
\n --blat=<file> blat rejected samples [default: NULL]
\n --pon=<file> pon from M2 [default: NULL]
\n --gnomad=<file> gnomad known variants [default: NULL]
\n --abra2=<file> result from running abra2 [default: NULL]" -> doc

library("docopt")
opts <- docopt(doc)

MergeVCFRecall=function(vcffile,outputstr, genome='hg19',runMode='Paired', blat=NULL, pon=NULL, gnomad=NULL, abra2=NULL){
  library(vcfR)
  library(matrixStats)
  print('Loading vcf sample of interest... ')
  BaseVcf=read.vcfR(vcffile)
  ponvcf=read.vcfR(pon)
  
  ## Set the location and index values?
  BaseVcfID=paste(BaseVcf@fix[ ,1], BaseVcf@fix[ ,2], BaseVcf@fix[ ,4])
  ponVcfID=paste(ponvcf@fix[ ,1], ponvcf@fix[ ,2], ponvcf@fix[ ,4])
  #head(blatRm)
  print('Process PoN data ...')
  head(ponvcf)
  # reannote the FILTER
  lx2=which(BaseVcfID %in% ponVcfID) # does not take into account indels
  BaseVcf@fix[lx2 ,7]=paste(BaseVcf@fix[lx2 ,7], "MutectPoN")
  rm(ponvcf)
  
  if (blat!="NULL"){
    print('Process blat samples....')
    blatRm=read.delim(blat, sep="\t")
    blatID=paste(blatRm$Chromosome, blatRm$Start_position, blatRm$Reference_Allele)
    lx1=which(BaseVcfID %in% blatID)
    BaseVcf@fix[lx1 ,7]=paste(BaseVcf@fix[lx1 ,7], "blat")
  }
  
  if (gnomad!="NULL"){
    print('Process gnomad data...')
    gnomadvcf=read.vcfR(gnomad)
    gnomadVcfID=paste(gnomadvcf@fix[ ,1], gnomadvcf@fix[ ,2], gnomadvcf@fix[ ,4])
    head(gnomadvcf)
    lx3=which(BaseVcfID %in% gnomadVcfID)
    BaseVcf@fix[lx3 ,7]=paste(BaseVcf@fix[lx3 ,7], "gnomad")
    rm(gnomadvcf)
  }
  
  if (abra2!="NULL"){
    print('Process abra2 data...')
    abra2vcf=read.vcfR(abra2)
    ##BaseVcfID2=paste(BaseVcf@fix[ ,1], BaseVcf@fix[ ,2])
    abra2VcfID=paste(abra2vcf@fix[ ,1], abra2vcf@fix[ ,2], abra2vcf@fix[ ,4])
    ##abra2VcfID2=paste(abra2vcf@fix[ ,1], abra2vcf@fix[ ,2])
    lx4=which(!BaseVcfID%in%abra2VcfID)
    BaseVcf@fix[lx4 ,7]=paste(BaseVcf@fix[lx4 ,7], "abra2Fail")
    ## Do a quick count of the number of columns which have an entry
    lx5=which(!abra2VcfID%in%BaseVcfID)
    tmp1=abra2vcf[lx5]
    tmp1@fix[ ,7]="abra2PASS"
    BaseVcf=rbind2(BaseVcf, tmp1)
    rm(abra2vcf)
  }
  sprintf('There are %s variants processed. Writing to file....', nrow(BaseVcf))
  #change the PASS annotations where it's missing
  BaseVcf@fix[ ,7]=gsub("PASS ", "",BaseVcf@fix[ ,7] )
  # write the output to file
  write.vcf(BaseVcf, file=paste(outputstr, ".judgement.vcf.gz", sep=""))
  # write the output to maf
  ## write the row max to file:
  tmp2=extract.gt(BaseVcf, element='AD', as.numeric=F)
  tmpB=sapply(strsplit(tmp2, ","), function(x) as.numeric(x[2]))
  tmpA=sapply(strsplit(tmp2, ","), function(x) as.numeric(x[1]))
  MatA=matrix(tmpA, ncol=ncol(tmp2), nrow=nrow(tmp2))
  MatB=matrix(tmpB, ncol=ncol(tmp2), nrow=nrow(tmp2))
  if (runMode=="Paired"){
    RefVal=rowMaxs(MatA[ , seq(2, ncol(MatA), by=2)])
    TVal=rowMaxs(MatB[ , seq(2, ncol(MatA), by=2)])
  }else{
    RefVal=rowMaxs(MatA)
    TVal=rowMaxs(MatB)    
  }
  MafFile=data.frame('Hugo_Symbol'=NA, 'Chromosome'=BaseVcf@fix[ ,1], 'Start_position'=BaseVcf@fix[ ,2],
                      't_ref_count'=RefVal, 't_alt_count'=TVal, 'Reference_Allele'=BaseVcf@fix[ ,4],
                      'Tumor_Seq_Allele2'=BaseVcf@fix[ ,5], 'judgement'='KEEP')
  write.table(MafFile, file=paste(outputstr, ".judgement.maf", sep=""), row.names = F, quote = F)
}

print(opts)

MergeVCFRecall(opts$vcffile, opts$outputstr, opts$genome,opts$runMode, opts$blat, opts$pon, opts$gnomad, opts$abra2)

