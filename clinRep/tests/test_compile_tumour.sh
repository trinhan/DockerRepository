#!/bin/bash
## Test compilation of the output report for Tumour
## Note if there is an error locally go into Rstudio and type Sys.getenv("RSTUDIO_PANDOC")
# Tumour report status
sampleName="NADOM22-001T"
ploidyTar="example_data/NADOM22-001T.optimalClusters.txt"

echo RUNNING TEST COMPILATION FOR ${sampleName} Tumour Mode

cp $ploidyTar .

# Functions to repopulate the .yaml file
export snvsummary="${sampleName}variantSummary.filt.maf"
export snvcancer="${sampleName}Cosmic.filt.maf"
export snvvus="${sampleName}VUS.filt.maf"
export snvdrug="${sampleName}Drug.filt.maf"
export snvhallmark="${sampleName}Pathway.filt.maf"
export snvacmg="${sampleName}ACMG.filt.maf"
export svSummary="${sampleName}.SV.SummaryTable.txt"
export svFull="${sampleName}.SV.full.filt.maf"
export svSplit="${sampleName}.SV.split.filt.maf"
export svACMG="${sampleName}.SV.acmg.filt.maf"
export cnvSummary="${sampleName}.CNV.SummaryTable.txt"
export cnvFull="${sampleName}.CNV.full.filt.maf"
export cnvSplit="${sampleName}.CNV.split.filt.maf"
export cnvACMG="${sampleName}.CNV.acmg.filt.maf"
export titanParams="${sampleName}.optimalClusters.txt"

rm -f final.yml temp.yml
( echo "cat > final.yml <<EOF"; cat example_data/${sampleName}.yaml; echo "EOF";)>temp.yml
. temp.yml

cp Template_Somatic_Report.Rmd ${sampleName}_Somatic_Report.Rmd

Path2="${sampleName}_Somatic_Report.Rmd"
echo run the compilation
Rscript -e "Sys.setenv(RSTUDIO_PANDOC='/Applications/RStudio.app/Contents/MacOS/quarto/bin');rmarkdown::render('./${sampleName}_Somatic_Report.Rmd')"   