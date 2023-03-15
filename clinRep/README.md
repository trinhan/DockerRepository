# README

Version 2.7

## Table of Contents

1. [Dockerfile](#dockerfile)  
2. [SNV data](#snv-data)
3. [SV data](#sv-data)
4. [CNV data](#cnv-data)
5. [Running Report](#clinical-report)
6. [Testing](#testing)  

<a name="headers"/>

## Dockerfile

The dockerfile can be built using:

```
docker build -t trinhanne/clin_report_annot:v2.7 .
```

## File Preparation

There are different series of scripts to prepare the inputs prior to embedding into the report (This does not run computation, it mainly displays the results)

### SNV data

SNV data is prepared using the following scripts. Note:

* The input can either be a vep annotated file, or a maf file following oncokb annotation. If the former, an extra function needs to be called
* The same process is run for both Germline or Tumour. The `--runMode` tag in `vepVCF2maf.R` will account for this 
* Smaller dependent files are included in the `annotFiles` folder if small enough, otherwise they are available on google cloud

#### Inputs

**Sample Dependent Inputs**

* vep annotated vcf.gz or oncokb .maf
* sampleName - to rename entries
* caseName - the actual ID of the sample in vcf or maf. Important to filter based on the tumour instead of normal

**Dependent Annotation Files**

File | Description | Download
-----|-------------|----------
AAlist | Amino acid 3 letter to single letter conversion for HGVSp_Short | annotFiles
AddList | Additional User defined list of genes to search for | annotFiles/AddList: Immunotherapy, DNADamage,chondrosarcoma
ACMG | A list of ACMG v3.0 genes | annotFiles, downloaded from PMID 34012068, Supp Table 1 
columnIDs | Parameters of interest to keep in final report and renaming | annotFiles
cosmicMut | List of consensus mutations from Cosmic | Download from cosmic, available on GCP
cosmicGenes | List of cosmic consensus genes | annotFiles or Download from cosmic 
MSigDB | MSigDB hallmark pathways to identify genes of interest, h.all.v7.4.symbols.gmt provided | annotFiles or download from [MsigDB](www.gsea-msigdb.org/gsea/msigdb/human/collections.jsp)
pathwayList | key search terms to return from MSigDB | annotFiles
pfam | protein domains | annotFiles, downloaded from [EBI](https://ftp.ebi.ac.uk/pub/databases/Pfam/releases/)
pirsf | protein domains | downloaded from [PIR](https://proteininformationresource.org/pirwww/download/ftpcenter.shtml) and grep function on ">"
scoringRubrik | scoring Matrix to define the 5 tiers | annotFiles

#### Execution

```
#Step 0: If input is a vcf file from VEP:
Rscript vepVCF2maf.R --vcffile ${vcffile} --outputfile ${outputmaf} --sampleName ${SampleName} --canonical T --runMode 'Germline' --AAlist annotFiles/AAlist.csv
#Step 1: Otherwise, If input is a maf file from Oncokb:
Rscript DBAnnotations.R --maffile ${outputmaf} --outputfile dba.maf --sampleName $SampleName --cosmicMut $cosmicMut --cosmicGenes $cosmicGenes --MSigDB $MsigDBAnnotation --pfam annotFiles/Pfam-A.clans.33.1.tsv --pirsf annotFiles/pirsfinfo_22072015.cat
paste vep.maf dba.maf > $Output
Rscript SummarizeVariants.R --maffile $Output --outputname $SampleName --caseName $SampleName --caddscore 10 --Ncallerthresh 1 --AddList $AddList --columnEntries annotFiles/ColumnIDs.csv
Rscript FilterVariants.R --maffile ${SampleName}variantsCoding.filt.maf --scoringRubrik annotFiles/scoringMatrix.csv --outputname $SampleName --ACMG annotFile/ACMG_v3.0.csv --pathwayList $pathwayTerms --pathwayFile annotFile/PathwayList.csv --gnomadcutoff 0.1 --onlyCoding T --pathogenic T
```

#### Outputs

* Scoring Matrix with criteria fulfilled for each gene
* Filtered gene lists:
	* all variants
	* ACMG only
	* Cosmic/Cancer
	* Pathway of interest
	* VUS
	* Drug Related


### SV data

* Input SV file is prepared the same way regardless of whether it was run on Tumour or Normal
* The assumption is the Tumour was run in single-sample mode anyway so it does not make a difference

#### Inputs

**Sample Dependent Inputs**

* AnnotSV annotated manta file
* sampleName - to rename entries
* caseName - the actual ID of the sample in vcf or maf. Important to filter based on the tumour instead of normal

**Dependent Annotation Files**

File | Description | Download
-----|-------------|----------
AddList | Additional User defined list of genes to search for | annotFiles/AddList: Immunotherapy, DNADamage,chondrosarcoma
GTex | Database GTex - to identify genes which are expressed in tissue of interest |
MSigDB | MSigDB hallmark pathways to identify genes of interest, h.all.v7.4.symbols.gmt provided | annotFiles or download from [MsigDB](www.gsea-msigdb.org/gsea/msigdb/human/collections.jsp)
pathwayList | key search terms to return from MSigDB | annotFiles

#### Execution

```
Rscript SummarizeAnnotSV.R --AnnotSVtsv ~{inputSV} --outputname ~{sampleName} --MSigDB ~{MsigDBAnnotation} \
        --GTex ~{GTex} --CosmicList ~{cosmicGenes} ~{"--AddList " + AddList} --pathwayList "$keyWd" --ACMGCutoff ~{ACMGcutoff} --Tissue "$TissueWd" --CNV ~{CNV} --PASSfilt TRUE
```

### CNV File preparation

## Clinical Report

## Testing


