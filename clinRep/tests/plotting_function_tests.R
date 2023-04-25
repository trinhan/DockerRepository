# Test required R functions

source("R/FilterSVs.R")
source("R/FilterCNVs.R")
#source("R/CreateGRangesData.R")

library(tools)
##############################
## Test 1: Filter SVs function
##############################
SVData="example_data/ATEST_NA12878_small_1x.SV.formated.tsv"
VAF=0.25
acmg=3
# Test 1: works as intended
Test1=tryCatch(expr={
                    FilterSVs(SVData,2, acmg); 
                    print("TEST1: FilterSV from germline checks PASSED")
                    }, 
                    error = function(e){
                      print("TEST1: FilterSV from germline FAILED")
                      })

##############################
## Test 2: Filter CNVs function from germline
##############################
CNVData="example_data/ATEST_NA12878_small_1x.gCNV.annotSV.tsv.gz"
CNlow=1
CNhigh=3
acmg=3
# Test 1: works as intended
Test1=tryCatch(expr={
  FilterCNVs(CNVData,CNlow,CNhigh, acmg); 
  print("TEST2: FilterCNV germline checks PASSED")
}, 
error = function(e){
  print("TEST2: FilterCNV germline FAILED")
})

##############################
## Test 3: Filter CNVs function from Funcotator
##############################
CNVData="example_data/ATEST_NA12878_small_1x.CNV.formated.tsv"
CNlow=1
CNhigh=3
acmg=3
# Test 1: works as intended
Test1=tryCatch(expr={
  FilterCNVs(CNVData,CNlow,CNhigh, acmg); 
  print("TEST2: FilterCNV germline checks PASSED")
}, 
error = function(e){
  print("TEST2: FilterCNV germline FAILED")
})