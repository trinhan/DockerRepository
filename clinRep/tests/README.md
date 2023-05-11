# README

Version 2.7

## List of Tests

Run the following tests in the order below to check everything works:

1. `bash tests/test_SNV_input.sh`
	- Test 1: ATEST_NA12878_small_1x.GRCh38_vep.vcf.gz | Germline | Vep annotated
	- Test 2: ER099_MEL4_Tempus_oncokb.maf.gz | Germline | Oncokb Annotated
	- Test 3: A_TEST_PAIRED.GRCh38_vep.vcf.gz | Tumour-Paired | Vep annotated
2. `bash tests/test_SV_input.sh`
	- Test 1: ER002_MPM1T_145_155bp.called.seg.funcotated.tsv | Tumour | GATK CNV
	- Test 2: ATEST_NA12878_small_1x.gCNV.annotSV.tsv.gz | Germline | GATK CNV 
	- Test 3: TEST_TUM_ER002_MPM1.Manta.annotSV.tsv.gz | Tumour | Manta
	- Test 4: ATEST_NA12878_small_1x.Manta.annotSV.tsv.gz | Germline | Manta
3. `Rscript test/plotting_functions_tests.R`
	- Test 1: CreateGRangesData
4. `bash test/test_compile_germline.sh`
	- Test 1: ATEST_NA12878_small_1x | Germline to generate the output file


