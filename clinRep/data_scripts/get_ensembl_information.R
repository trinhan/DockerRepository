## Do this for entire transcriptome/genome and match?
## download from https://www.genenames.org/download/custom/, 29/07/2022
allTx=read.delim("~/Downloads/gene_names1.txt", sep="\t")
ensembl <- useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
genesAll <- toGRanges(getBM(attributes=c('chromosome_name', 'start_position', 'end_position', 'hgnc_symbol'),
                         filters = 'hgnc_symbol', values =allTx$Approved.symbol, mart = ensembl))
save(genesAll, file="~/gitLibs/DockerRepository/clinRep/annotFiles/ensembl_granges_genes_290722.RData")
