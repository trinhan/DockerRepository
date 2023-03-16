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

## SNV data

SNV data is prepared using the following scripts. Note:

* The input can either be a vep annotated file, or a maf file following oncokb annotation. If the former, an extra function needs to be called
* The same process is run for both Germline or Tumour. The `--runMode` tag in `vepVCF2maf.R` will account for this 
* Smaller dependent files are included in the `annotFiles` folder if small enough, otherwise they are available on google cloud

#### Inputs

**Sample Dependent Inputs**

* vep annotated vcf.gz or oncokb .maf
* sampleName - to rename entries
* caseName - the actual ID of the sample in vcf or maf. Important to filter based on the tumour instead of normal
* pathwayList - which pathway of interest to narrow-in on (e.g. Immunotherapy, DNADamage)


**Dependent Annotation Files**

File | Description | Download
-----|-------------|----------
AAlist | Amino acid 3 letter to single letter conversion for HGVSp_Short | annotFiles
AddList | Additional User defined list of genes to search for | annotFiles/AddList: Immunotherapy, DNADamage,chondrosarcoma
ACMG | A list of ACMG v3.0 genes | annotFiles, downloaded from PMID 34012068, Supp Table 1 
columnIDs | Parameters of interest to keep in final report and renaming | annotFiles
cosmicMut | List of Cancer Mutation Census Data from Cosmic | Download from [cosmic](https://cancer.sanger.ac.uk/cosmic/download), available on GCP
cosmicGenes | List of cosmic consensus genes | annotFiles or Download from [cosmic Gene Census](https://cancer.sanger.ac.uk/census)
MSigDB | MSigDB hallmark pathways to identify genes of interest, h.all.v7.4.symbols.gmt provided | annotFiles or download from [MsigDB](www.gsea-msigdb.org/gsea/msigdb/human/collections.jsp)
pfam | protein domains | annotFiles, downloaded from [EBI](https://ftp.ebi.ac.uk/pub/databases/Pfam/releases/)
pirsf | protein domains | downloaded from [PIR](https://proteininformationresource.org/pirwww/download/ftpcenter.shtml) and grep function on ">"
scoringRubrik | scoring Matrix to define the 5 SNV tiers | annotFiles
pathwayFile | key search terms to return from MSigDB | annotFiles 

**Filtering Thresholds**

Below is just a suggestion of what to use. Note that thresholds for caddscore, gnomad may be applicable to only some of the outputs - see the scoring Rubrik to check which "SNV Tier" considers it

Parameter | Summary |Threshold or Value 
---------------|--------------|-------------
canonical |  If multiple transcripts are present in vep output, select the canonical one  | TRUE 
caddscore | Predicted deleteriousness by CADD. Higher score is more deleterious | 10 as threshold of top 10\% deleteriousness
Ncallerthresh | Minimum number of callers needed to support a variant | Suggestion of 2 or higher
gnomadcutoff | Population frequency cut-off (if used in scoring criteria) | 0.1 - expect most of the variants to be synergistic with the treatment of interest and not deleterious in its own right.
onlyCoding | Limit the outputs to only variants in coding regions and UTRs | T
pathogenic | Prediction of pathogenicity based on either cadd, clinvar, IMPACT, polyphen | TRUE - limits the output lists

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


## SV data

* Input SV file is prepared the same way regardless of whether it was run on Tumour or Normal
* The assumption is the Tumour was run in single-sample mode anyway so it does not make a difference

#### Inputs

**Sample Dependent Inputs**

* AnnotSV annotated manta file
* sampleName - to rename entries
* caseName - the actual ID of the sample in vcf or maf. Important to filter based on the tumour instead of normal
* Tissue - the tissue in which the tumour is found 
* pathwayTerm - which pathway of interest to narrow-in on (e.g. Immunotherapy, DNADamage)
* germline - whether the input is germline or not

**Dependent Annotation Files**

File | Description | Download
-----|-------------|----------
AddList | Additional User defined list of genes to search for | annotFiles/AddList: Immunotherapy, DNADamage, chondrosarcoma
cosmicGenes | List of cosmic consensus genes | annotFiles or Download from [cosmic Gene Census](https://cancer.sanger.ac.uk/census)
GTex | Database GTex - to identify genes which are expressed in tissue of interest |
MSigDB | MSigDB hallmark pathways to identify genes of interest, h.all.v7.4.symbols.gmt provided | annotFiles or download from [MsigDB](www.gsea-msigdb.org/gsea/msigdb/human/collections.jsp)
pathwayList | key search terms to return from MSigDB | annotFiles

**Filtering Thresholds**

Below is just a suggestion of what to use. Note that thresholds for caddscore, gnomad may be applicable to only some of the outputs - see the scoring Rubrik to check which "SNV Tier" considers it

Parameter | Summary |Threshold or Value 
---------------|--------------|-------------
ACMGCutoff |  AnnotSV ACMG scoring of the variant | 3. Note that most high-scoring SVs (5) are whole chromosomal deletions 
SRfilter | Minimum number of spanning reads to support a variant | 5? depends on the sequencing depth
PRfilter | Minimum number of paired reads to support a variant | 2 
PASSfilt | Retain only variants which are considered "PASS" by manta if this information is present | TRUE

#### Execution

```
Rscript SummarizeAnnotSV.R --AnnotSVtsv ~{inputSV} --outputname ~{sampleName} --germline T --MSigDB ${MsigDBAnnotation} \
        --GTex ${GTex} --CosmicList ${cosmicGenes} --AddList ${AddList} --pathwayTerm ${pathwayTerm} --pathwayList ${pathwayList} --ACMGCutoff ~{ACMGcutoff} --Tissue ${Tissue} --CNV FALSE --SRfilter 5 --PRfilter 2 --PASSfilt TRUE
```

## CNV Data

Note that CNVs were run through different pipelines if run in germline mode vs in tumour mode

* Germline Mode: `GATK` followed by ` AnnotSV`
* Tumour Mode: GATK followed by Funcotator

As a result, there are two different preparation steps currently. TDB whether to fuse these approaches

### Germline Mode

The process is identical to the SV pre-processing. Note that the only difference is a "CNV" vs "SV" tag within the SummarizeAnnotSV.R function, and SR and PR filter values do not need to be specified 

#### Execution

```
Rscript SummarizeAnnotSV.R --AnnotSVtsv ~{inputSV} --outputname ~{sampleName} --germline T --MSigDB ${MsigDBAnnotation} \
        --GTex ${GTex} --CosmicList ${cosmicGenes} --AddList ${AddList} --pathwayTerm ${pathwayTerm} --pathwayList ${pathwayList} --ACMGCutoff ${ACMGcutoff} --Tissue ${Tissue} --CNV T --PASSfilt TRUE
```

### Tumour Mode

**Inputs**

* A tsv file from funcotator
* The annotation files are the same as for Germline

NOTE: The main difference here is that the input file is different. Review the code to see if we can combine this all into one script

#### Execution

```
Rscript AnnotateTumCNV.R --tsv ${inputCNV} --outputname ${sampleName} --MSigDB ${MsigDBAnnotation} \
        --GTex ${GTex} --CosmicList ${cosmicGenes} --AddList  ${AddList} --pathwayTerm ${pathwayTerm} --pathwayList ${pathwayList} --Tissue ${Tissue}
```

## Clinical Report

## Testing


