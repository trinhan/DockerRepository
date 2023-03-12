#!/usr/bin/Rscript
"This script tidies an input maffile by collapsing genotypes and filtering out variants most likely to have a consequence

TierE. ACMG - Genes within the ACMG guidelines
TierA. COSMIC - Known oncogenes
TierB. PATHWAY - Relevant pathways to the expected disease mechanism, and/or treatment resistance mechanisms
TierC. DRUG - Genes with clinvar indications of 
TierD. VUS - All other genes with a functional consequence

usage: \n SummarizeVariants.R [--maffile=<file> --outputname=<string> --caseName=<string> --caddscore=<float> --gnomadcutoff=<float> --Ncallerthresh=<int> --AddList=<file> --columnEntries=<file>]
\n options:
\n --maffile=<file> maffile that has been annotated using vep & oncokb.
\n --outputname=<string> output string [default: output.tsv]
\n --caseName=<string> name of the case samplex
\n --caddscore=<float> [default: 20] score for cadd cut-off
\n --gnomadcutoff=<float> [default: 0.1]
\n --Ncallerthresh=<int> [default: 0] minimum number of callers required to output the variant
\n --AddList=<file> Additional user defined genes [default=NULL]
\n --columnEntries=<file> Import a table of columns to export and renamed output " -> doc

library("docopt", quietly = T)
opts <- docopt(doc)

############################################
#1. Load in all required libraries silently
############################################

suppressMessages(library(data.table, quietly = T))
suppressMessages(library(dplyr, quietly = T))
writeLines('Running SummarizeVariants \n Read in original File.. \n')

#####################################################################
#2. Read in first line of maffile and extract colnames wtih on the total reads(nTot), genotype (GT), VAF and nAlt reads for each caller
#####################################################################  

mydnames=fread(opts$maffile, sep="\t", nrows=0)
gnames=colnames(mydnames)[grep("nTot", colnames(mydnames))]
gnames2=colnames(mydnames)[grep("^GT", colnames(mydnames))]
gnames3=colnames(mydnames)[grep("VAF", colnames(mydnames))]
gnames4=colnames(mydnames)[grep("nAlt", colnames(mydnames))]

#####################################################################
#3. Read in file of column headers to extract
#####################################################################  

ColNames=read.csv(opts$columnEntries, header=F)
## check which colnames matches the original file
midx=match(ColNames$V1, colnames(mydnames))
mNames=colnames(mydnames)[na.omit(midx)]

#####################################################################
#4. Read in maffile, only with column names of interest (in step 2 and 3)
#####################################################################  

InputData=fread(opts$maffile, sep="\t", select=c(mNames, gnames, gnames2, gnames3, gnames4))

######################################################################
#5. Remove HLA genes (hypervariable) and sep annotations from genotypes
#######################################################################
rm1=c(grep("HLA", InputData$Hugo_Symbol))
AnnotOnly=InputData[setdiff(1:nrow(InputData), rm1), 1:length(mNames) ]
ValsOnly=InputData[setdiff(1:nrow(InputData), rm1), (length(mNames)+1):ncol(InputData) ]
newNames=ColNames$V2[which(!is.na(midx))]
colnames(AnnotOnly)=ColNames$V2[which(!is.na(midx))]

#####################################################
##6. rename columns if they don't exist - if bypass onocokb, some inputs may not exist
###################################################
if ( !"ClinVar.Disease"%in%newNames & "ClinVar_Trait"%in%newNames){
    print('rename ClinVar Traits as ClinVar Disease')
    colnames(AnnotOnly)[which(newNames=="ClinVar_Trait")]="ClinVar.Disease"
  }
  
if (!"ClinVar.Sig"%in%newNames & "ClinSig"%in%newNames){
     print('rename column ClinSig as ClinVar.Sig')
     colnames(AnnotOnly)[which(newNames=="ClinSig")]="ClinVar.Sig"
  }
  
  if (!"gnomAD.AF"%in%newNames & "gnomad_MAX_AF"%in%newNames){
    print('rename column gnomad_MAX_AF as gnomAD.AF')
    colnames(AnnotOnly)[which(newNames=="gnomad_MAX_AF")]="gnomAD.AF"
  }

#####################################################
##7. reset the CADD scoring threshold if column does not exist
###################################################

  if (!"CADD"%in%newNames){
    print('reset CADD cut-off, not part of selection criteria')
    opts$caddscore=0
  }

#####################################################
##8. Consider max value of gnomad and MGRB and assign in $AF_max. Set unknown AF to 0
###################################################
  
  if (!"MGRB.AF"%in%newNames & "gnomAD.AF"%in%newNames){
    print('gnomad only')
    AnnotOnly$AF_max=AnnotOnly$gnomAD.AF
    AnnotOnly$AF_max[which(is.na(AnnotOnly$AF_max))]=0
  } else if ("MGRB.AF"%in%newNames & "gnomAD.AF"%in%newNames){
    print('both gnomad and mgrb exist, use the higher AF')
    AnnotOnly$gnomAD.AF[which(is.na(AnnotOnly$gnomAD.AF))]=0
    AnnotOnly$MGRB.AF[which(is.na(AnnotOnly$MGRB.AF))]=0
    AnnotOnly$AF_max=ifelse(AnnotOnly$gnomAD.AF>AnnotOnly$MGRB.AF, AnnotOnly$gnomAD.AF, AnnotOnly$MGRB.AF)
  } else if ("MGRB.AF"%in%newNames & !"gnomAD.AF"%in%newNames){
    print('MGRB only')
    AnnotOnly$AF_max=AnnotOnly$MGRB.AF
    AnnotOnly$AF_max[which(is.na(AnnotOnly$AF_max))]=0
  }

#####################################################
##9. Clean up ClinVar Annotations and Protein 
###################################################

print('Parse ClinVar and Protein annotations')
AnnotOnly$ClinVar.Disease=gsub("not_specified", "", AnnotOnly$ClinVar.Disease)
AnnotOnly$ClinVar.Disease=gsub("not_provided", "", AnnotOnly$ClinVar.Disease)
AnnotOnly$ClinVar.Disease=gsub("&", " ", AnnotOnly$ClinVar.Disease)
AnnotOnly$ProteinDomain=ifelse(is.na(AnnotOnly$pfam), AnnotOnly$pirsf, AnnotOnly$pfam)

#####################################################
##10. Modify disease frequency summary
###################################################

print('Modify disease frequency text')
#e.g. upper_aerodigestive_tract/carcinoma/squamous_cell_carcinoma=20/935=2.14%;central_nervous_system/glioma/astrocytoma_Grade_I=10/12=83.33%;soft_tissue/haemangioblastoma/NS=22/25=88%"
#to "upper_aerodigestive_tract SCC=2.14%;central_nervous_system glioma astrocytoma_Grade_I=83.33%;soft_tissue haemangioblastoma=88%"

#remove "NS" - number of samples in text
AnnotOnly$Cancer.Disease.Freq=gsub("\\/NS", "", AnnotOnly$Cancer.Disease.Freq)
#remove "/carcinoma/" as non-specific and usually accompanied by SCC or NSCC
AnnotOnly$Cancer.Disease.Freq=gsub("\\/carcinoma\\/", " ", AnnotOnly$Cancer.Disease.Freq)
#abbreviate squamous cell carcinoma
AnnotOnly$Cancer.Disease.Freq=gsub("squamous_cell_carcinoma", "SCC", AnnotOnly$Cancer.Disease.Freq)
#abbreviate non small cell carcinoma
AnnotOnly$Cancer.Disease.Freq=gsub("non_small_cell_carcinoma", "NSCC", AnnotOnly$Cancer.Disease.Freq)
#remove the number of cases, and only consider percentages
AnnotOnly$Cancer.Disease.Freq=gsub("=[0-9]*\\/[0-9]*", "", AnnotOnly$Cancer.Disease.Freq)
#remnove slashes and replace with space
AnnotOnly$Cancer.Disease.Freq=gsub("\\/", " ", AnnotOnly$Cancer.Disease.Freq)


#####################################################
##11. Compute summary stats on number of callers, average depth etc
###################################################

print('Extract genotypes and callers specific to the sample')
# Extract nTotal
Nx=ValsOnly %>% select(all_of(gnames))
Nx=data.matrix(data.frame(Nx))
Nxb=Nx[ , grep(opts$caseName, colnames(Nx))]
  
## Extract GT
GT=ValsOnly %>% select(all_of(gnames2))
GT=data.frame(GT)
GTb=GT[, grep(opts$caseName, colnames(GT))]

## Extract VAF
VAF=ValsOnly %>% select(all_of(gnames3))
VAF=data.frame(VAF)
VAFb=VAF[ , grep(opts$caseName, colnames(VAF)) ]

## Extract nAlt
nAlt=ValsOnly %>% select(all_of(gnames4))
nAlt=data.frame(nAlt)
nAltb=nAlt[ , grep(opts$caseName, colnames(nAlt)) ]

##If a sample has nAlt=0, substitute as NA
tx=which(nAltb==0, arr.ind = T)
nAltb[tx]=NA
VAFb[tx]=NA
GTb[tx]=NA

## put in a catch for single callers
  if (is.null(ncol(Nxb))){
    print('note only 1 caller used')
    # header name - caller.nTot.Name
    HeaderV=strsplit(colnames(Nx),"\\.")
    HeaderV=sapply(HeaderV, function(x) x[3])
    # set number of callers to 1
    AnnotOnly$Ncall=1
    # list the name of the caller
    AnnotOnly$Ncallers=unlist(HeaderV)[1]
    # Set the depth, Genotype and VAF 
    AnnotOnly$Depth=Nxb
    AnnotOnly$Genotype=GTb
    AnnotOnly$VAF=VAFb
  }else{
    HeaderV=strsplit(colnames(Nxb),"\\.")
    HeaderV=sapply(HeaderV, function(x) x[3])
    # Calculate the av depth from all callers
    AvDepth=ceiling(rowMeans(Nxb, na.rm=T))
    # count the number of callers for each row
    Ncall=rowSums(sign(Nxb), na.rm = T)
    # obtain the names of the callers which do not return NA, and separate with ,
    Nx2=sapply(1:ncol(Nxb), function(x) ifelse(!is.na(Nxb[,x ]), HeaderV[x], NA))
    Ncallers=sapply(1:nrow(Nx2), function(x) paste(na.omit(Nx2[x, ]),collapse=","))
    # assign values
    AnnotOnly$Ncallers=Ncallers
    AnnotOnly$Ncall=Ncall
    AnnotOnly$Depth=AvDepth
    # compute the VAF and round to 3 decimal places
    VAF2=rowMeans(VAFb, na.rm = T)
    AnnotOnly$VAF=round(VAF2, 3)
    # Collapse the Genotypes - this assumes they will all be similar?
    # Take the most common genotype for each row, exclude NA values
    GTcons2=sapply(1:nrow(GTb), function(x) names(which.max(table(t(GTb[x, ])))))
    AnnotOnly$Genotype=GTcons2
  }

# filter out variants if they dont meet the minimum caller requirement
sprintf('Keep only variants called by more than %s callers', opts$Ncallerthresh)
AnnotOnly=AnnotOnly[which(AnnotOnly$Ncall>opts$Ncallerthresh), ]

#####################################################
##12. Add annotation if gene is in user list
###################################################

  AnnotOnly$UserGeneList=NA
  if (!is.null(opts$AddList)){
    # assume the file doesn't have a header
    UserTab=read.csv(opts$AddList, header=F)
    # extract only the first line as a vector
    UserTab=UserTab[ ,1]
    # indicate it is in the gene list with a value 1
    int2=which(AnnotOnly$SYMBOL%in%UserTab)
    AnnotOnly$UserGeneList[int2]=1
  }
  
#####################################################
## 13. Add extra annotation on predicted pathogenicicty:
## * coding consequence * pathogenicity filter based on clinvar, polyphen,IMPACT and/or CADD 
###################################################

  ## Annotate whether the variant could have a functional consequence
  ConsequenceVals=c("missense", "nonsense", "frameshift", "splice", "UTR", "inframe", "stop", "NMD")
  Nxgrep=sapply(ConsequenceVals, function(x) grep(x, AnnotOnly$Consequence))
  AnnotOnly$ConsB=NA
  AnnotOnly$ConsB[unlist(Nxgrep)]=1
  
  # Add a potential filter based on pathogenicity
  AnnotOnly$Pathogenicity=NA
  # select based on clinVar
  inClinVar=unique(unlist(sapply(c("pathogenic", "risk_factor", "drug_response", "protective"), function(x) grep(x, AnnotOnly$ClinVar.Sig, ignore.case=T))))
  # select based on polyphen
  inPolyPhen=grep("damaging", AnnotOnly$PolyPhen)
  # select in "IMPACT"
  inIMPACT=which(AnnotOnly$IMPACT=="HIGH")
  # select based on CADD
  if (as.numeric(opts$caddscore)>0){
    inCADD=which(AnnotOnly$CADD> as.numeric(opts$caddscore))
    idx=unique(c(inClinVar, inPolyPhen, inIMPACT,inCADD))
  }else{
    idx=unique(c(inClinVar, inPolyPhen, inIMPACT))
  }
  AnnotOnly$Pathogenicity[idx]=1

#####################################################
## 14. Filter an output table for the next steps
## save to file: all Variants, filtered and summary table
###################################################

  print('Generate filteredlist based on gnomad, protein consequence and Pathogenicity')
  nidx=which(AnnotOnly$AF_max<=as.numeric(opts$gnomadcutoff) & AnnotOnly$ConsB==1 & AnnotOnly$Pathogenicity==1)
  FilteredTable=AnnotOnly[nidx,  ]
  sprintf("Summarize variants with gnomad %s Consequence coding, Expanded list contains %s variants", as.numeric(opts$gnomadcutoff), length(nidx))
    
  ## Summary Stats
  DFValues=data.frame(Nvar=nrow(AnnotOnly))
  
  print('Finished! Writing out files')
  FilteredTable <- apply(FilteredTable,2,as.character)
  AnnotOnly <- apply(AnnotOnly,2,as.character)
  write.table(FilteredTable, file=paste(opts$outputname, "variantsCoding.filt.maf", sep=""), sep = "\t", row.names = F,  quote = F)
  write.table(DFValues, file=paste(opts$outputname, "variantSummary.filt.maf", sep=""), sep = "\t", row.names = F,  quote = F)
  write.table(AnnotOnly, file=paste(opts$outputname, "variantsAll.maf", sep=""), sep = "\t", row.names = F,  quote = F)
 