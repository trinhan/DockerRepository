# Oncokb

## Summary

This container uses Oncokb (current version 3.2.4) as downloaded from this github repository: https://github.com/oncokb/oncokb-annotator

It contains additional custom Rscripts for the following functions:

* **vepVCF2maf.R** - This script converts an output vcf.gz from VEP into a maf file. Note that this function is identical to that listed in the clinReport folder

## Example Usage (within the Docker container)

For 1 sample within 1 vcf file:

```
docker run -v ${pwd}:/mnt -it trinhanne/oncokb:v3.2.4Renv
Rscript opt/vepVCF2maf.R --vcffile ${vcffile} --outputfile ${Output} --sampleName ${SampleName} --canonical T --runMode 'Germline' --AAlist annotFiles/AAlist.csv
python3 /oncokb/MafAnnotator.py -i ${Output} -o "${SampleName}_oncokb.maf" -c "clinannot.txt" -b ${token} -q HGVSp_Short
```

Note that you will need register with oncokb https://www.oncokb.org/account/register in order to obtain an access token ${token} 

## Example Usage (within WDL)

Please see the workflow: https://github.com/trinhan/WDLPipelines/blob/main/oncokb.wdl

Note the corresponding json https://github.com/trinhan/WDLPipelines/blob/main/json/oncokb_inputs.json is outdated 

