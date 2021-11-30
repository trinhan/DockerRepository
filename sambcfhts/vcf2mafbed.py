'''
Script for saving vcf file to maf for blat
Also for outputting a bed file
'''
### run python strelka_allelic_count_snv.py ${STRELKA2_SNV_PASS} $STRELKA2_SNV_REFORMATTED_VCF "${caseName}" "${ctrlName}" "Method"


import allel
import pandas as pd
import numpy as np
from argparse import ArgumentParser


def _process_vcf_file(invcf, outMaf, runMode, outBed, bpPad):
	''' Function to read the file and output maf and bed file
		by adding alelic depth (AD) and alelic fraction (AF) to it
	'''
	vcf=allel.read_vcf(invcf, fields=['samples', 'variants/CHROM','variants/POS', 'variants/REF', 'variants/ALT', 'calldata/AD'])
	tmpVal=vcf['calldata/AD']
	if runMode=="Paired":
		lsp=np.arange(start=1, stop=tmpVal.shape[1], step=2)
	elif runMode=="TumOnly":
		lsp=np.arange(start=1, stop=tmpVal.shape[1], step=1)
	tmpV2=tmpVal[:,lsp, 0]
	Mref=tmpV2.max(axis=1)
	tmpVA=tmpVal[:,lsp, 1]
	Malt=tmpVA.max(axis=1)
	newdata=pd.DataFrame({'Hugo_Symbol':["na"]*len(Mref), 'Chromosome':vcf['variants/CHROM'], 'Start_position':vcf['variants/POS'],'t_ref_count':Mref, 't_alt_count':Malt, 'Reference_Allele':vcf['variants/REF'], 'Tumor_Seq_Allele2':vcf['variants/ALT'][:, 0],'judgement':["KEEP"]*len(Mref)})
	posBuffer=vcf['variants/POS']-bpPad
	posBuffer[(posBuffer<0)]=0
	posBuffer=[str(int) for int in posBuffer]
	posBufferUp=vcf['variants/POS']+bpPad
	posBufferUp=[str(int) for int in posBufferUp]
	varray=np.vstack((vcf['variants/CHROM'], posBuffer, posBufferUp))
	t2=['\t'.join(varray[:,x]) for x in range(1, len(Mref))]
	textfile = open(outBed, "w")
	for element in t2:
		textfile.write(element + "\n")
	textfile.close()
	newdata.to_csv(outMaf, sep='\t', index=False)

def main():
# The main argument parser
	parser = ArgumentParser(description="Export vcf entries as maf (for blat) and bed file (for realignment) ")
	parser.add_argument('IN_VCF', help='name of the input vcf file')
	parser.add_argument('OUT_MAF',help='name of the output maf file')
	parser.add_argument('OUT_BED',help="name of the output bed file")
	parser.add_argument('bpPad', help="integer of padding (e.g. 100 bases) up or downstream for bed file", type=int, default=150)
	parser.add_argument('runMode', help="was the variant calculated in paired (Paired) or tumor only (TumOnly)? ", type=str, default="Paired")
	args = parser.parse_args()

	if(args):
		_process_vcf_file(args.IN_VCF, args.OUT_MAF, args.OUT_BED, args.bpPad)
	else:
		parser.print_help()


if __name__ == "__main__":
	main()




