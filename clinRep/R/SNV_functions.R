## SNV functions

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
  colnames(adv)=paste("nAlt", colnames(adv), sep="-")
  colnames(dp)=paste("nTot", colnames(dp), sep="-")
  colnames(vaf)=paste("VAF", colnames(vaf), sep="-")
  colnames(gt)=paste("GT", colnames(gt), sep="-")
  
  VAFd=cbind(adv, dp, vaf, gt)
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