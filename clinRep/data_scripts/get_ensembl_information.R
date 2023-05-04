## Do this for entire transcriptome/genome and match?
## download from https://www.genenames.org/download/custom/, 29/07/2022
allTx=read.delim("./gene_names_290722.txt", sep="\t")
ensembl <- useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
genesAll <- toGRanges(getBM(attributes=c('chromosome_name', 'start_position', 'end_position', 'hgnc_symbol'),
                         filters = 'hgnc_symbol', values =allTx$Approved.symbol, mart = ensembl))
seqlevelsStyle(genesAll)<-"UCSC"
## filter to retain only the genes within chr1:22 and X,Y,MT
chrNames=c(paste("chr", c(1:22), sep=""), "chrX", "chrY", "chrMT")
genesAll=genesAll[which(seqnames(genesAll)%in%chrNames)]
save(genesAll, file="~/gitLibs/DockerRepository/clinRep/annotFiles/ensembl_granges_genes_290722_chr1_22_XYMT.RData")
