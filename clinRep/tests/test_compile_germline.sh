#!/bin/bash
## Test compilation of the output report for Germline
## Note if there is an error locally go into Rstudio and type Sys.getenv("RSTUDIO_PANDOC")
# Germline report status
sampleName="ATEST_NA12878_small_1x"
ploidyTar="example_data/ATEST_NA12878_small_1x.case-contig-ploidy-calls.tar.gz"

echo RUNNING TEST COMPILATION FOR ${sampleName} Germline Mode

mkdir ploidytmp
tar xf ${ploidyTar} -C ploidytmp
mv ploidytmp/*/contig_ploidy.tsv .

rm -rf ploidytmp

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

rm -f final.yml temp.yml
( echo "cat > final.yml <<EOF"; cat example_data/${sampleName}.yaml; echo "EOF";)>temp.yml
. temp.yml

cp Template_Germline_Report.Rmd ${sampleName}_Germline_Report.Rmd

Path2="${sampleName}_Germline_Report.Rmd"
echo run the compilation
Rscript -e "Sys.setenv(RSTUDIO_PANDOC='/Applications/RStudio.app/Contents/MacOS/quarto/bin');rmarkdown::render('./${sampleName}_Germline_Report.Rmd')"   

#Rscript -e 'library(rmarkdown);Sys.setenv(RSTUDIO_PANDOC="/Applications/RStudio.app/Contents/MacOS/quarto/bin"); rmarkdown::render(eval("templateRmd/${sampleName}_Germline_Report.Rmd"))'   
