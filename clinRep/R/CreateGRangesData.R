## Inputs:
## 1. ASV- full table: require CHROM, START, END
## 2. Col - colour scale to use
## 3. SVmode - True or F
## 3. Ncheck - If in SV mode, this value to merge together variants within this distance

CreateGRangesData=function(ASV, Col, SVmode, Ncheck=500){
  if (SVmode){
    #####################################
    ## set up the Granges object for SV
    #####################################   
    gMat <- GRanges(seqnames=paste0("chr",ASV[ ,"CHROM"]),ranges=IRanges(start=ASV[, "START"], end=ASV[, "END"]),names=ASV$SV_ID, 
                    path=ASV$PathogenicSource, AF=ASV$MeanVAF)
    # search distances between all starts and select smaples for merge
    tx1=as.matrix(dist(start(gMat)))
    tx2=as.matrix(dist(end(gMat)))
    ReviewThese=which(tx1<Ncheck & tx2<Ncheck, arr.ind=T)
    ReviewThese=ReviewThese[which(ReviewThese[ ,1]<ReviewThese[ ,2]), ]
    CheckChr=which(seqnames(gMat)[ReviewThese[ ,1]]==seqnames(gMat)[ReviewThese[ ,2]])
    if (length(CheckChr)>0 ){
      #tx3=ReviewThese[ which(ReviewThese[ ,1]!=ReviewThese[ ,2]), ] ## Is this actually needed?
      tx4=ReviewThese[CheckChr, ]
      for (i in 1:nrow(tx4)){
        ## Merge the samples
        end(gMat[tx4[i, 1]])=end(gMat[tx4[i, 2]]) 
        gMat$AF[tx4[i, 1]]=mean(gMat$AF[tx4[i, ]])
      }
      gMat=gMat[-tx4[,2 ]]
    }
  } else {
    #############################################
    # set up the Granges object for CNV purpose
    #############################################
    gMat <- GRanges(seqnames=paste0("chr",ASV[, "CHROM"]), ranges=IRanges(start=ASV[, "START"], end=ASV[, "END"]), 
                       path=ASV$PathogenicSource, AF=ifelse(is.na(ASV$Pop.Freq), 1, 1-ASV$Pop.Freq), gainloss=ifelse(ASV$CN>2, "gain", "loss"))
  }
  
  gMat$AFcol=brewer.pal(6, Col)[cut(gMat$AF, c(0,0.25, 0.5, 0.8, 0.9, 1))]
  
  return(gMat)  
}