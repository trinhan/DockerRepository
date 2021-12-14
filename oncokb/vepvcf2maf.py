'''
Script for saving vcf file to maf for blat
Also for outputting a bed file
'''
### run python strelka_allelic_count_snv.py ${STRELKA2_SNV_PASS} $STRELKA2_SNV_REFORMATTED_VCF "${caseName}" "${ctrlName}" "Method"


import allel
import pandas as pd
import numpy as np
from argparse import ArgumentParser


def _process_vcf_file(invcf, vepvcf, outMaf, outBed, bpPad, runMode):
	''' Function to read the file and output maf and bed file
		by adding alelic depth (AD) and alelic fraction (AF) to it
	'''
	vcf=allel.read_vcf(invcf, fields=['samples', 'variants/CHROM','variants/POS', 'variants/REF', 'variants/ALT', 'variants/FILTER_PASS','calldata/AD'])
	tmpVal=vcf['calldata/AD']
	tmpV2=tmpVal[:,:, 0]
	tmpVA=tmpVal[:,:, 1]
	tx1=tmpVal.shape[1]
	V2_names=[x + "_Ref_Counts" for x in vcf['samples']]
	VA_names=[x + "_Alt_Counts" for x in vcf['samples']]
	t2=pd.DataFrame(data=tmpV2, columns=V2_names)
	ta=pd.DataFrame(data=tmpVA, columns=VA_names)
	newdata=pd.DataFrame({'Chromosome':vcf['variants/CHROM'], 'Start_position':vcf['variants/POS'], 'Reference_Allele':vcf['variants/REF'], 'Tumor_Seq_Allele2':vcf['variants/ALT'][:, 0], 'Filter':vcf['variants/FILTER_PASS'] })
	df_c = pd.concat([newdata.reset_index(drop=True), t2, ta], axis=1)
	print(df_c)
	## Now read the vep annotation
	vep=allel.read_vcf(vepvcf)
	print(vep)
	newdata.to_csv(df_c, sep='\t', index=False)

def main():
# The main argument parser
	parser = ArgumentParser(description="Export vcf entries as maf (for blat) and bed file (for realignment) ")
	parser.add_argument('IN_VCF', help='name of the input vcf file')
	parser.add_argument('VEP_VCF', help='name of the input vcf file')
	parser.add_argument('OUT_MAF',help='name of the output maf file')
	parser.add_argument('OUT_BED',help="name of the output bed file")
	parser.add_argument('bpPad', help="integer of padding (e.g. 100 bases) up or downstream for bed file", type=int, default=150)
	parser.add_argument('runMode', help="was the variant calculated in paired (Paired) or tumor only (TumOnly)? ", type=str, default="Paired")
	args = parser.parse_args()

	if(args):
		_process_vcf_file(args.IN_VCF, args.VEP_VCF, args.OUT_MAF, args.OUT_BED, args.bpPad, args.runMode)
	else:
		parser.print_help()


if __name__ == "__main__":
	main()




