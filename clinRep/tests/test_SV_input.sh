#!/bin/bash
## Run SV/CNV pre-processing for different scenarios

#######################
# Dependencies
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
pathwayTerm="DNADamage"
pathwayList="annotFiles/PathwayList.csv"
columnEntries="annotFiles/SV_column_IDs.csv"
Tissue="skin"

#############################
# Run Specific Cases
# -- toggle below to select which tests to run
# 0 is off/false; 1 is on/true
#############################
TEST1=1
TEST2=0
TEST3=0
TEST4=0
TEST5=0
TEST6=1

#############################
# Test 1: CNV Funcotator
#############################

if [ $TEST1 -eq 1 ]; then
  echo "RUNNING TEST 1 - TUMOUR CNV"
  inputSV="example_data/ER002_MPM1T_145_155bp.called.seg.funcotated.tsv"
  sampleName="TEST_TUM_ER002_MPM1"
  CNV=TRUE
  germline=FALSE
  Funco=TRUE

  Rscript scripts/SummarizeAnnotSV.R --tsv ${inputSV} --outputname ${sampleName} --germline ${germline} --PASSfilt FALSE --MSigDB ${MsigDBAnnotation} --GTex ${GTex} --CosmicList ${cosmicGenes} --CNV ${CNV} --AddList ${AddList} --pathwayTerm ${pathwayTerm} --pathwayList ${pathwayList} --Tissue ${Tissue} --funco ${Funco}
  Rscript R/FilterCNVsFunco.R --tsv ${sampleName}.CNV.formated.tsv --outputname ${sampleName} --CNlow 0 --CNhigh 4

  RESULT=$?
  if [ $RESULT -eq 0 ]; then
    echo TEST1 PASSED
  else
    echo TEST1 FAILED
  fi
fi

#############################
# Test 2: CNV Germline
#############################
if [ $TEST2 -eq 1 ]; then
inputSV="example_data/ATEST_NA12878_small_1x.gCNV.annotSV.tsv.gz"
sampleName="ATEST_NA12878_small_1x"
CNV=TRUE
germline=TRUE
PassFilt=FALSE
Funco=FALSE

echo "RUNNING TEST 2 - GERMLINE CNV"
Rscript scripts/SummarizeAnnotSV.R --tsv ${inputSV} --outputname ${sampleName} --germline ${germline} --PASSfilt ${PassFilt} --MSigDB ${MsigDBAnnotation} --GTex ${GTex} --CosmicList ${cosmicGenes} --CNV ${CNV} --AddList ${AddList} --pathwayTerm ${pathwayTerm} --pathwayList ${pathwayList} --Tissue ${Tissue} --ACMGCutoff 2 --funco ${Funco}
Rscript R/FilterSVs.R --tsv ${sampleName}.CNV.formated.tsv --outputname ${sampleName} --mode CNV --CNlow 0 --CNhigh 4 --ColsIDs ${columnEntries}
RESULT=$?
if [ $RESULT -eq 0 ]; then
  echo TEST2 PASSED
else
  echo TEST2 FAILED
fi
fi

#############################
# Test 3: SV call - tumour
#############################
if [ $TEST3 -eq 1 ]; then
inputSV="example_data/NADOM22-001T.Manta.annotSV.tsv.gz"
sampleName="NADOM22-001T"
CNV=FALSE
germline=FALSE
PassFilt=TRUE
Funco=FALSE

echo "RUNNING TEST 3 - SV TUMOUR"
Rscript scripts/SummarizeAnnotSV.R --tsv ${inputSV} --outputname ${sampleName} --germline ${germline} --PASSfilt ${PassFilt} --MSigDB ${MsigDBAnnotation} --GTex ${GTex} --CosmicList ${cosmicGenes} --CNV ${CNV} --AddList ${AddList} --pathwayTerm ${pathwayTerm} --pathwayList ${pathwayList} --Tissue ${Tissue} --ACMGCutoff 2 --SRfilter 3 --PRfilter 0 --funco ${Funco} &&
Rscript R/FilterSVs.R --tsv ${sampleName}.SV.formated.tsv --outputname ${sampleName} --mode SV --VAF 0.2 --ACMGcutoff 5 --ColsIDs ${columnEntries}
RESULT=$?
if [ $RESULT -eq 0 ]; then
  echo TEST3 PASSED
else
  echo TEST3 FAILED
fi
fi

#############################
# Test 4: SV call - germline
#############################
if [ $TEST4 -eq 1 ]; then
inputSV="example_data/ATEST_NA12878_small_1x.Manta.annotSV.tsv.gz"
sampleName="ATEST_NA12878_small_1x"
CNV=FALSE
germline=TRUE
PassFilt=FALSE
Funco=FALSE

echo "RUNNING TEST 4 - SV GERMLINE"
Rscript scripts/SummarizeAnnotSV.R --tsv ${inputSV} --outputname ${sampleName} --germline ${germline} --PASSfilt ${PassFilt} --MSigDB ${MsigDBAnnotation} --GTex ${GTex} --CosmicList ${cosmicGenes} --CNV ${CNV} --AddList ${AddList} --pathwayTerm ${pathwayTerm} --pathwayList ${pathwayList} --Tissue ${Tissue} --ACMGCutoff 2 --SRfilter 3 --PRfilter 0 --funco ${Funco} && 
Rscript R/FilterSVs.R --tsv ${sampleName}.SV.formated.tsv --outputname ${sampleName} --mode SV --VAF 0.2 --ACMGcutoff 5
RESULT=$?
if [ $RESULT -eq 0 ]; then
  echo TEST4 PASSED
else
  echo TEST4 FAILED
fi
fi

#############################
# Test 5: GATK CNV Tumour + AnnotSV
#############################

if [ $TEST5 -eq 1 ]; then
  echo "RUNNING TEST 5 - TUMOUR CNV on AnnotSV"
  inputSV="example_data/TEST_TUM_ER002_MPM1.gatkSomaticCNV.annotSV.tsv.gz"
  sampleName="TEST_TUM_ER002_MPM1"
  CNV=TRUE
  germline=FALSE
  Funco=FALSE

  Rscript scripts/SummarizeAnnotSV.R --tsv ${inputSV} --outputname ${sampleName} --germline ${germline} --PASSfilt FALSE --MSigDB ${MsigDBAnnotation} --GTex ${GTex} --CosmicList ${cosmicGenes} --CNV ${CNV} --AddList ${AddList} --pathwayTerm ${pathwayTerm} --pathwayList ${pathwayList} --Tissue ${Tissue}  --ACMGCutoff 5 --funco ${Funco}&&
  Rscript R/FilterSVs.R --tsv ${sampleName}.CNV.formated.tsv --outputname ${sampleName} --mode CNV --CNlow -1 --CNhigh 1

  #Rscript R/FilterCNVsFunco.R --tsv ${sampleName}.CNV.formated.tsv --outputname ${sampleName} --CNlow -1 --CNhigh 1

  RESULT=$?
  if [ $RESULT -eq 0 ]; then
    echo TEST5 PASSED
  else
    echo TEST5 FAILED
  fi
fi

#############################
# Test 6: GATK CNV Tumour + AnnotSV
#############################

if [ $TEST6 -eq 1 ]; then
  echo "RUNNING TEST 6 - TUMOUR CNV on AnnotSV - real data"
  inputSV="example_data/NAMDOM22-002T.gatkSomaticCNV.annotSV.tsv.gz"
  sampleName="NAMDOM22-002T"
  CNV=TRUE
  germline=FALSE
  Funco=FALSE

  Rscript scripts/SummarizeAnnotSV.R --tsv ${inputSV} --outputname ${sampleName} --germline ${germline} --PASSfilt FALSE --MSigDB ${MsigDBAnnotation} --GTex ${GTex} --CosmicList ${cosmicGenes} --CNV ${CNV} --AddList ${AddList} --pathwayTerm ${pathwayTerm} --pathwayList ${pathwayList} --Tissue ${Tissue}  --ACMGCutoff 5 --funco ${Funco} &&
  Rscript R/FilterSVs.R --tsv ${sampleName}.CNV.formated.tsv --outputname ${sampleName} --mode CNV --CNlow -1 --CNhigh 1

  #Rscript R/FilterCNVsFunco.R --tsv ${sampleName}.CNV.formated.tsv --outputname ${sampleName} --CNlow -1 --CNhigh 1

  RESULT=$?
  if [ $RESULT -eq 0 ]; then
    echo TEST6 PASSED
  else
    echo TEST6 FAILED
  fi
fi