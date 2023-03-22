#!/usr/bin/Rscript
"Usage: \n 
  DBAnnotations.R (--maffile=<file> --cosmicMut=<file> --cosmicGenes=<file> --MSigDB=<file> --pfam=<file> --pirsf=<file> ) [--outputfile=<string> --sampleName=<string> ]
  
options:
 --maffile=<file>
 --outputfile=<string> output file name [default: output.txt]
 --sampleName=<string> If sample list has more than 1 sample, include a csv file with all identifiers, otherwise set as 'NULL' [default: NULL]
 --cosmicMut=<file> File containing list of Cosmic Mutations
 --cosmicGenes=<file> File of Cosmic Cancer Genes
 --MSigDB=<file> File containing MsigDB gene set information 
 --pfam=<file> File containing pfam annotations [default: NULL]
 --pirsf=<file> File containing psirf annotations [default: NULL] " -> doc

library("docopt", quietly = T)
opts <- docopt(doc, help=TRUE, version='1.0.0')

############################################
#1. Load in all required libraries silently
############################################
suppressMessages(library(data.table, quietly = T))
suppressMessages(library(GSEABase, quietly = T))
suppressMessages(library(dplyr, quietly = T))
  
#####################################################################
#2. Read in the relevant headers in the output maf file to save space
#####################################################################  
AllData=fread(opts$maffile, sep="\t", stringsAsFactors = F, select=c("Hugo_Symbol","CHROM","POS", "REF", "ALT", "HGVSp_Short","HGVSc", "DOMAINS"))
sprintf('Running DBAnnotations...\n reading in maffile...\n Table dimensions %s rows by %s col', nrow(AllData), ncol(AllData))
  
#########################################
#3. Annotate with the pathways in MSigDB
#########################################
## a. read in the MSigSB file and replace the "HALLMARK" at the beginning of each gene set
print('annotate with MSigDB..')
PathInH=getGmt(con=opts$MSigDB, geneIdType=SymbolIdentifier(), collectionType=BroadCollection(category="h"))
Mx2=geneIds(PathInH)
Ids2=names(PathInH)
Ids2=gsub("HALLMARK_", "", Ids2)
names(Mx2)=Ids2
Mx3=stack(Mx2)
## b. Find all the unique genes in our maffile
GS=unique(AllData$Hugo_Symbol)
## c. Find all the pathways which feature the gene 
Tx2=sapply(GS, function(x) paste(Mx3$ind[which(Mx3$values==x)], collapse=" "))
## d. Match and assign the pathway list to the original DF 
Nm2=match(AllData$Hugo_Symbol, GS)
AllData$HallmarkPathways=Tx2[Nm2]
## e. Finished statement, and remove the MSigDB item to save space 
sprintf('%s Genes are now annotated with MSigDB pathway', length(Nm2))
rm(PathInH)

#########################################
#4. Annotate with Cosmic Mutational List
#########################################
  
print('annotate with Cosmic')
## a. read in required columns from the Mutation List
CosmicD=fread(opts$cosmicMut, sep="\t", select=c("GENE_NAME","POS","Mutation AA","ONC_TSG","CGC_TIER","DISEASE", "CLINVAR_TRAIT","MUTATION_SIGNIFICANCE_TIER"))
# b. Match based on chromosomal position  
aax1=paste(AllData$Hugo_Symbol, AllData$POS)
bbx1=paste(CosmicD$GENE_NAME, CosmicD$POS)
m1=match(aax1, bbx1)
# b. Match based on mutation AA
ax1=paste(AllData$Hugo_Symbol, AllData$HGVSp_Short)
bx1=paste(CosmicD$GENE_NAME, CosmicD$`Mutation AA`)
mmx1=match(ax1, bx1)
m1[which(!is.na(mmx1))]=mmx1[which(!is.na(mmx1))]
# c. Filter out the required columns from Cosmic and append to original Data Frame
CosmicD=CosmicD[m1,c("Mutation AA", "POS","ONC_TSG","CGC_TIER","DISEASE", "CLINVAR_TRAIT","MUTATION_SIGNIFICANCE_TIER") ]
colnames(CosmicD)=paste("CMC", colnames(CosmicD), sep=".")
AllData=cbind(AllData, CosmicD)
# d. remove the Cosmic data.frame to save space
rm(CosmicD)

############################################
#5. Annotate with Consensus Cosmic Gene List
############################################
# a. load the file and match with our existing gene list
print('check Cosmic Genes')
Cosmic2=read.csv(opts$cosmicGenes)
tx2=match(AllData$Hugo_Symbol, Cosmic2$Gene.Symbol)
# b. Add the following information: Tier
AllData$CMC.Cancer_Gene_Tier=NA
AllData$CMC.Cancer_Gene_Tier[which(!is.na(tx2))]=Cosmic2$Tier[na.omit(tx2)]
# c. If the gene is a hallmark gene, add this information in the Tier column 
tx3=match(AllData$Hugo_Symbol, Cosmic2$Gene.Symbol[which(Cosmic2$Hallmark=="Yes")])
AllData$CMC.Cancer_Gene_Tier[which(!is.na(tx3))]=paste("Hallmark", AllData$CMC.Cancer_Gene_Tier[which(!is.na(tx3))])

#########################################
#6. Annotate with Consensus Cosmic Gene List
#########################################

print('Add protein domain annotations from pfam and psird')
# a. read in data
pfamI=read.delim(opts$pfam, header=F, sep="\t")
pirsfinfo=read.delim(opts$pirsf, sep=")", header=F, quote="")
# b. pirsf data frame: add a new column which takes input after (
pirsfinfo$V3=NA
pirsfinfo$V3=substr(pirsfinfo$V1, 2, 12)
# c. Match the "PF00000" values in the domains column with the pfam dataset
print('add pfam')
n2=regmatches( AllData$DOMAINS,regexec("PF[0-9]{5}", AllData$DOMAINS))
n2pfam=unlist({n2[sapply(n2, length)==0] <- NA; n2})
pfamO=as.character(pfamI$V5[match(n2pfam, pfamI$V1)])
# c. Match the "PIRSF000000" values in the domains column with the pirsf dataset
print('add pirsf')
n2=regmatches( AllData$DOMAINS, regexec("PIRSF[0-9]{6}", AllData$DOMAINS))
n2=unlist({n2[sapply(n2, length)==0] <- NA; n2})
pirsfO=as.character(pirsfinfo$V2[match(n2,pirsfinfo$V3)])
  
#########################################
# 7. Save data and write to file
#########################################
# Note the output will be row concatenated to the existing file - can omit the first 4 lines of information as it will be repeated
AllData=cbind(AllData, pfamid=n2pfam, pfam=pfamO, pirsf=pirsfO, pirsfid=n2)
sprintf('write to file. Table dimensions %s row by %s columns', nrow(AllData), ncol(AllData))
write.table(AllData[ ,-c(1:4)], file=opts$outputfile, sep="\t", quote=F,row.names = F)

print('DBANNOTATIONS COMPLETE')
