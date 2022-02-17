#!/usr/bin/Rscript
"usage: \n annotateProteins.R [--maffile=<file> --outputfile=<string> --pfam=<file> --pirsf=<file> ]
\n options:\n --maffile=<file> maffile with protein domain annotation.
\n --outputfile=<string> output directory
\n --pfam=<file> File containing pfam annotations [default: NULL]
\n --pirsf=<file> File containing psirf annotations [default: NULL] " -> doc

library("docopt")
opts <- docopt(doc)
#opts

## vcf to maf format
## following annotation with VEP104
## run with the options:
## --domains --af_gnomad --check_existing --fields CANONICAL,Consequence,Codons,Amino_acids,Gene,SYMBOL,Feature,EXON,PolyPhen,SIFT,Protein_position,BIOTYPE,Feature_type,cDNA_position,CDS_position,Existing_variation,DISTANCE,STRAND,CLIN_SIG,LoF_flags,LoF_filter,LoF,RadialSVM_score,RadialSVM_pred,LR_score,LR_pred,CADD_raw,CADD_phred,Reliability_index,DOMAINS,HGVSc,HGVSp,gnomAD_AF,PHENO

annotateProteins=function(maffile,outputfile,pfam, pirsf){
  AllData=read.delim(maffile, sep="\t")
  print('converting proteins')
  ## pfam entry
  pfamI=read.delim(pfam, header=F, sep="\t")
  pirsfinfo=read.delim(pirsf, sep=")", header=F, quote="")
  pirsfinfo$V3=substr(pirsfinfo$V1, 2, 12)
  ##pirsfinfo$V1=gsub(">", "", pirsfinfo$V1)
  
  n2=regmatches( AllData[ ,match("DOMAINS", colnames(AllData))], 
                 regexec("PF[0-9]{5}", AllData[ ,match("DOMAINS", colnames(AllData))]))
  n2pfam=unlist({n2[sapply(n2, length)==0] <- NA; n2})
  pfamO=as.character(pfamI$V5[match(n2pfam, pfamI$V1)])
  
  n2=regmatches( AllData[ ,match("DOMAINS", colnames(AllData))], 
                 regexec("PIRSF[0-9]{6}", AllData[ ,match("DOMAINS", colnames(AllData))]))
  n2=unlist({n2[sapply(n2, length)==0] <- NA; n2})
  pirsfO=as.character(pirsfinfo$V2[match(n2,pirsfinfo$V3)])
  
  AllData=cbind(AllData, pfamid=n2pfam, pfam=pfamO, pirsf=pirsfO, pirsfid=n2)
  write.table(AllData, file=outputfile, sep="\t", row.names=F, quote=F)
}

annotateProteins(opts$maffile, opts$outputfile, opts$pfam, opts$pirsf)
