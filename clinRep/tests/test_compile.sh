#!/bin/bash
## Test compilation of the output report for Germline
## Note if there is an error locally go into Rstudio and type Sys.getenv("RSTUDIO_PANDOC")
# Germline report status
sampleName="ATEST_NA12878_small_1x"

echo RUNNING TEST COMPILATION FOR ${sampleName} Germline Mode

# Functions to repopulate the .yaml file
export snvsummary="${sampleName}variantSummary.filt.maf"
export snvcancer="${sampleName}Cosmic.filt.maf"
export snvvus="${sampleName}VUS.filt.maf"
export snvdrug="${sampleName}Drug.filt.maf"
export snvhallmark="${sampleName}Pathway.filt.maf"
export snvacmg="${sampleName}ACMG.filt.maf"
export cnvAnnot="${sampleName}.CNV.formated.tsv"
export svAnnot="${sampleName}.SV.formated.tsv"
export SVsplit="${sampleName}.SV.split.formated.tsv"

rm -f final.yml temp.yml
( echo "cat > final.yml <<EOF"; cat example_data/${sampleName}.yaml; echo "EOF";)>temp.yml
. temp.yml

cp Template_Germline_Report.Rmd ${sampleName}_Germline_Report.Rmd

Path2="${sampleName}_Germline_Report.Rmd"
echo run the compilation
Rscript -e "Sys.setenv(RSTUDIO_PANDOC='/Applications/RStudio.app/Contents/MacOS/quarto/bin');rmarkdown::render('./${sampleName}_Germline_Report.Rmd')"   

#Rscript -e 'library(rmarkdown);Sys.setenv(RSTUDIO_PANDOC="/Applications/RStudio.app/Contents/MacOS/quarto/bin"); rmarkdown::render(eval("templateRmd/${sampleName}_Germline_Report.Rmd"))'   
