#!/bin/bash
## Test all the tests in this script

Sample="/Users/anntri/Desktop/test_space_clin_rep/test_MEL4_tempus/ER099_MEL4.GRCh38_vep.vcf.gz"
SampleName="1x_test"
Output="1x_test_output.maf"
ACMG="/Users/anntri/Documents/ER_pilot/annotations/ACMG_v3.0.csv"
AAlist="/Users/anntri/gitLibs/DockerRepository/clinRep/annotFiles/AminoAcid_table.csv"
cosmicMut="/Users/anntri/Downloads/CosmicMut_export_v92_pos_annot.tsv"
cosmicGenes="/Users/anntri/Documents/ER_pilot/annotations/CosmicGeneCensus_all03142021.csv"
MsigDBAnnotation="/Users/anntri/Documents/ER_pilot/annotations/h.all.v7.4.symbols.gmt"
pfam="/Users/anntri/Documents/ER_pilot/annotations/Pfam-A.clans.33.1.tsv"
pirsf="/Users/anntri/Documents/ER_pilot/annotations/pirsfinfo_22072015.cat"
GTex="/Users/anntri/Documents/ER_pilot/annotations/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz"
AddList="/Users/anntri/Documents/ER_pilot/annotations/DNADamage_50_genes.csv"
ScoreMatrix="/Users/anntri/gitLibs/DockerRepository/clinRep/scoringMatrix.csv"
pathwayTerms="Immunotherapy"
pathwayList="/Users/anntri/gitLibs/DockerRepository/clinRep/annotFiles/PathwayList.csv"
columnEntries="/Users/anntri/gitLibs/DockerRepository/clinRep/annotFiles/ColumnIDs.csv"

Rscript vepVCF2maf.R --vcffile $Sample --outputfile vep.maf --sampleName $SampleName --canonical T --runMode Tumour --AAlist $AAlist
Rscript DBAnnotations.R --maffile vep.maf --outputfile dba.maf --sampleName $SampleName --cosmicMut $cosmicMut --cosmicGenes $cosmicGenes --MSigDB $MsigDBAnnotation --pfam $pfam --pirsf $pirsf 
paste vep.maf dba.maf > $Output
Rscript SummarizeVariants.R --maffile $Output --outputname $SampleName --caseName "ER099_MEL4" --caddscore 10 --Ncallerthresh 1 --AddList $AddList --columnEntries $columnEntries
Rscript FilterVariants.R --maffile ${SampleName}variantsCoding.filt.maf --scoringRubrik $ScoreMatrix --outputname $SampleName --ACMG $ACMG --pathwayList $pathwayTerms --pathwayFile $pathwayList --gnomadcutoff 0.1 --onlyCoding T --pathogenic T
