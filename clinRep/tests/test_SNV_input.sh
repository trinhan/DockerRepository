#!/bin/bash
## Run SNV processing tests for 3 different scenarios

#######################
# Other required inputs
#######################
cosmicMut="/Users/anntri/Downloads/CosmicMut_export_v92_pos_annot.tsv"
GTex="/Users/anntri/Documents/ER_pilot/annotations/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz"
ACMG="annotFiles/ACMG_v3.0.csv"
AAlist="annotFiles/AminoAcid_table.csv"
cosmicGenes="annotFiles/CosmicGeneCensus_all03142021.csv"
MsigDBAnnotation="annotFiles/h.all.v7.4.symbols.gmt"
pfam="annotFiles/Pfam-A.clans.33.1.tsv"
pirsf="annotFiles/pirsfinfo_22072015.cat"
AddList="annotFiles/AddList/DNADamage_50_genes.csv"
ScoreMatrix="annotFiles/scoringMatrix.csv"
pathwayTerms="DNADamage"
pathwayList="annotFiles/PathwayList.csv"
columnEntries="annotFiles/ColumnIDs.csv"

#############################
# Run Specific Cases
# -- toggle below to select which tests to run
# 0 is off/false; 1 is on/true
#############################
TEST1=0
TEST2=0
TEST3=0
TEST4=1
TEST5=1
TEST6=1

#######################
## Sample Options here
#######################
# 1. Germline WGS vep annotated - works
#####################
Sample="example_data/ATEST_NA12878_small_1x.GRCh38_vep.vcf.gz"
SampleName="ATEST_NA12878_small_1x"
Output="ATEST_NA12878_small_1x.maf"
runMode="Germline"

if [ $TEST1 -eq 1 ]; then
  echo "RUNNING TEST 1 - GERMLINE WGS VEP"
  Rscript scripts/vepVCF2maf.R --vcffile $Sample --outputfile vep.maf --sampleName $SampleName --canonical T --runMode $runMode --AAlist $AAlist &&
  Rscript scripts/DBAnnotations.R --maffile vep.maf --outputfile dba.maf --sampleName $SampleName --cosmicMut $cosmicMut --cosmicGenes $cosmicGenes --MSigDB $MsigDBAnnotation --pfam $pfam --pirsf $pirsf &&
  paste vep.maf dba.maf > $Output &&
  Rscript scripts/SummarizeVariants.R --maffile $Output --outputname $SampleName --caseName $SampleName --caddscore 10 --Ncallerthresh 1 --AddList $AddList --columnEntries $columnEntries &&
  Rscript scripts/FilterVariants.R --maffile ${SampleName}variantsCoding.filt.maf --scoringRubrik $ScoreMatrix --outputname $SampleName --ACMG $ACMG --pathwayList $pathwayTerms --pathwayFile $pathwayList --gnomadcutoff 0.1 --onlyCoding T --pathogenic T
  RESULT=$?
  rm vep.maf dba.maf

  if [ $RESULT -eq 0 ]; then
    echo TEST1 PASSED
  else
    echo TEST1 FAILED
  fi
fi

################################
# 2. Germline WGS oncokb annotated
###############################
Sample="example_data/ER099_MEL4_Tempus_oncokb.maf.gz"
SampleName="ER099_MEL4_Tempus"
Output="ER099_MEL4_germline.maf"
runMode="Germline"
gunzip -c $Sample > vep.maf

if [ $TEST2 -eq 1 ]; then
echo "Running Test 2 - GERMLINE WGS ONCOKB"
Rscript scripts/DBAnnotations.R --maffile vep.maf --outputfile dba.maf --sampleName $SampleName --cosmicMut $cosmicMut --cosmicGenes $cosmicGenes --MSigDB $MsigDBAnnotation --pfam $pfam --pirsf $pirsf &&
paste vep.maf dba.maf > $Output &&
Rscript scripts/SummarizeVariants.R --maffile $Output --outputname $SampleName --caseName $SampleName --caddscore 10 --Ncallerthresh 1 --AddList $AddList --columnEntries $columnEntries &&
Rscript scripts/FilterVariants.R --maffile ${SampleName}variantsCoding.filt.maf --scoringRubrik $ScoreMatrix --outputname $SampleName --ACMG $ACMG --pathwayList $pathwayTerms --pathwayFile $pathwayList --gnomadcutoff 0.1 --onlyCoding T --pathogenic T
RESULT=$?
   rm vep.maf dba.maf
if [ $RESULT -eq 0 ]; then
  echo TEST2 PASSED
else
  echo TEST2 FAILED
fi
fi

################################
# 3. Tumour Paired - works
# Note there is now an issue with reading like 1052
################################
 Sample="example_data/A_TEST_PAIRED.GRCh38_vep.vcf.gz"
 SampleName="TEST_TUM_ER002_MPM1"
 Output="A_TEST_PAIRED.maf"
 runMode="Tumour"

if [ $TEST3 -eq 1 ]; then 
echo "RUNNING TEST 3 - TUMOUR PAIRED VEP"
Rscript scripts/vepVCF2maf.R --vcffile $Sample --outputfile a_test_vep.maf --sampleName $SampleName --canonical T --runMode $runMode --AAlist $AAlist &&
Rscript scripts/DBAnnotations.R --maffile vep.maf --outputfile dba.maf --sampleName $SampleName --cosmicMut $cosmicMut --cosmicGenes $cosmicGenes --MSigDB $MsigDBAnnotation --pfam $pfam --pirsf $pirsf &&
paste a_test_vep.maf dba.maf > $Output &&
Rscript scripts/SummarizeVariants.R --maffile $Output --outputname $SampleName --caseName $SampleName --caddscore 10 --Ncallerthresh 1 --AddList $AddList --columnEntries $columnEntries &&
Rscript scripts/FilterVariants.R --maffile ${SampleName}variantsCoding.filt.maf --scoringRubrik $ScoreMatrix --outputname $SampleName --ACMG $ACMG --pathwayList $pathwayTerms --pathwayFile $pathwayList --gnomadcutoff 0.1 --onlyCoding T --pathogenic T
RESULT=$?
 rm a_test_vep.maf dba.maf
if [ $RESULT -eq 0 ]; then
  echo TEST3 PASSED
else
  echo TEST3 FAILED
fi
fi

################################
# 4. Tumour Only - works
################################
Sample="example_data/NADOM22-001.GRCh38_vep.vcf.gz"
SampleName="NADOM22-001T"
Output="NADOM22-001_Tum_Only.maf"
runMode="Tumour"

if [ $TEST4 -eq 1 ]; then 
echo "RUNNING TEST 4 - TUMOUR SINGLE SAMPLE VEP"
Rscript scripts/vepVCF2maf.R --vcffile $Sample --outputfile vep2.maf --sampleName $SampleName --canonical T --runMode $runMode --AAlist $AAlist &&
Rscript scripts/DBAnnotations.R --maffile vep2.maf --outputfile dba.maf --sampleName $SampleName --cosmicMut $cosmicMut --cosmicGenes $cosmicGenes --MSigDB $MsigDBAnnotation --pfam $pfam --pirsf $pirsf 
paste vep2.maf dba.maf > $Output &&
Rscript scripts/SummarizeVariants.R --maffile $Output --outputname $SampleName --caseName $SampleName --caddscore 10 --Ncallerthresh 1 --AddList $AddList --columnEntries $columnEntries &&
Rscript scripts/FilterVariants.R --maffile ${SampleName}variantsCoding.filt.maf --scoringRubrik $ScoreMatrix --outputname $SampleName --ACMG $ACMG --pathwayList $pathwayTerms --pathwayFile $pathwayList --gnomadcutoff 0.1 --onlyCoding T --pathogenic T
RESULT=$?
  rm vep2.maf dba.maf
if [ $RESULT -eq 0 ]; then
  echo TEST4 PASSED
else
  echo TEST4 FAILED
fi
fi


################################
# 4. Tumour Only - works
################################
Sample="example_data/NADOM22-003.GRCh38_vep.vcf.gz"
SampleName="NADOM22-003T"
Output="NADOM22-003_Tum_Only.maf"
runMode="Tumour"

if [ $TEST5 -eq 1 ]; then 
echo "RUNNING TEST 4 - TUMOUR SINGLE SAMPLE VEP"
Rscript scripts/vepVCF2maf.R --vcffile $Sample --outputfile vep2.maf --sampleName $SampleName --canonical T --runMode $runMode --AAlist $AAlist &&
Rscript scripts/DBAnnotations.R --maffile vep2.maf --outputfile dba.maf --sampleName $SampleName --cosmicMut $cosmicMut --cosmicGenes $cosmicGenes --MSigDB $MsigDBAnnotation --pfam $pfam --pirsf $pirsf 
paste vep2.maf dba.maf > $Output &&
Rscript scripts/SummarizeVariants.R --maffile $Output --outputname $SampleName --caseName $SampleName --caddscore 10 --Ncallerthresh 1 --AddList $AddList --columnEntries $columnEntries &&
Rscript scripts/FilterVariants.R --maffile ${SampleName}variantsCoding.filt.maf --scoringRubrik $ScoreMatrix --outputname $SampleName --ACMG $ACMG --pathwayList $pathwayTerms --pathwayFile $pathwayList --gnomadcutoff 0.1 --onlyCoding T --pathogenic T
RESULT=$?
  rm vep2.maf dba.maf
if [ $RESULT -eq 0 ]; then
  echo TEST4 PASSED
else
  echo TEST4 FAILED
fi
fi

################################
# 4. Tumour Only - works
################################
Sample="example_data/NAMDOM22-002.GRCh38_vep.vcf.gz"
SampleName="NAMDOM22-002T"
Output="NAMDOM22-002_Tum_Only.maf"
runMode="Tumour"

if [ $TEST6 -eq 1 ]; then 
echo "RUNNING TEST 4 - TUMOUR SINGLE SAMPLE VEP"
Rscript scripts/vepVCF2maf.R --vcffile $Sample --outputfile vep2.maf --sampleName $SampleName --canonical T --runMode $runMode --AAlist $AAlist &&
Rscript scripts/DBAnnotations.R --maffile vep2.maf --outputfile dba.maf --sampleName $SampleName --cosmicMut $cosmicMut --cosmicGenes $cosmicGenes --MSigDB $MsigDBAnnotation --pfam $pfam --pirsf $pirsf 
paste vep2.maf dba.maf > $Output &&
Rscript scripts/SummarizeVariants.R --maffile $Output --outputname $SampleName --caseName $SampleName --caddscore 10 --Ncallerthresh 1 --AddList $AddList --columnEntries $columnEntries &&
Rscript scripts/FilterVariants.R --maffile ${SampleName}variantsCoding.filt.maf --scoringRubrik $ScoreMatrix --outputname $SampleName --ACMG $ACMG --pathwayList $pathwayTerms --pathwayFile $pathwayList --gnomadcutoff 0.1 --onlyCoding T --pathogenic T
RESULT=$?
  rm vep2.maf dba.maf
if [ $RESULT -eq 0 ]; then
  echo TEST4 PASSED
else
  echo TEST4 FAILED
fi
fi