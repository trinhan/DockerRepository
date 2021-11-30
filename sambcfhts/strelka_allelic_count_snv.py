'''
Script for adding allelic depth and allelic fraction to the Strelka VCF file
Also for changing TUMOR and NORMAL to tumor and normal sample names
'''
### run python strelka_allelic_count_snv.py ${STRELKA2_SNV_PASS} $STRELKA2_SNV_REFORMATTED_VCF "${caseName}" "${ctrlName}" "Method"


from argparse import ArgumentParser

def _get_allelic_depth(fields, refbase, altbase, values):
	#d = dict(zip(fields.split(':'), values.split(':')))
	d=values.rstrip('\n').split(':')
	base_info_loc = {'A':4, 'C':5, 'G':6, 'T':7}
	
	ref_count = int(d[base_info_loc[refbase] ].split(',')[0])
	if (len(altbase)==1):	
		alt_count = int(d[base_info_loc[altbase] ].split(',')[0])
	else:
		alt_count=int(d[0])-ref_count
	
	return ref_count, alt_count 


def _compute_allelic_fraction(ref, alt):
	try:
		# Somatic allele frequency is $tier1AltCounts / ($tier1AltCounts + $tier1RefCounts)
		return float(alt) / (float(ref) + float(alt))
	except:
		print ('REF count (%s) or ALT count (%s) could not be converted into floats or both equal to 0.' %(ref, alt))
		return 0.0


def _process_vcf_file(inVcf,outVcf, TUMOR_NAME, NORMAL_NAME, VariantCaller):
	''' Function modifies Strelka indels VCF file 
		by adding alelic depth (AD) and alelic fraction (AF) to it
	'''
	tumor_idx = 10
	normal_idx = 9
	format_idx = 8
	info_idx = 7
	ref_idx=3
	alt_idx=4
	add_ad_af_description = True
	with open (outVcf, 'w') as writer:
		with open(inVcf, 'r') as reader:
			for line in reader:
				if line.startswith('##'):
					if line.startswith('##FORMAT') and add_ad_af_description:
						writer.write('##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">' + '\n')
						writer.write('##FORMAT=<ID=AF,Number=A,Type=Float,Description="Allele fractions of alternate alleles in the tumor">' + '\n')
						writer.write('##FORMAT=<ID=SBS,Number=A,Type=Float,Description="Strelka2 SNVSB score">' + '\n')
						writer.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">' + '\n')
						add_ad_af_description = False
					writer.write(line)
				elif line.startswith('#CHROM'):
					pieces = line.strip('\n').split('\t')
					#print pieces
					tmr_idx = pieces.index('TUMOR')
					nrml_idx = pieces.index('NORMAL')
					info_idx=pieces.index('INFO')
					pieces[tmr_idx] = TUMOR_NAME
					pieces[nrml_idx] = NORMAL_NAME
					writer.write('\t'.join(pieces) + '\n')
				else:
					pieces = line.strip('\n').split('\t')
					INFO = pieces[info_idx]
					if VariantCaller!="None":
						pieces[info_idx] =INFO + ';ME=' + str(VariantCaller)
					FORMAT = pieces[format_idx]
					NORMAL = pieces[normal_idx]
					TUMOR = pieces[tumor_idx]
					REFBASE = pieces[ref_idx]
					ALTBASE = pieces[alt_idx]
					NORMAL_AD_REF, NORMAL_AD_ALT = _get_allelic_depth(FORMAT, REFBASE, ALTBASE, NORMAL)
					NORMAL_AF = _compute_allelic_fraction(NORMAL_AD_REF, NORMAL_AD_ALT)
					TUMOR_AD_REF, TUMOR_AD_ALT = _get_allelic_depth(FORMAT, REFBASE, ALTBASE, TUMOR)
					TUMOR_AF = _compute_allelic_fraction(TUMOR_AD_REF, TUMOR_AD_ALT)
					SNVSV=INFO.strip('\n').split(';')
					snsvs2=SNVSV[11]
					snsvs3=snsvs2.strip('\n').split('=')
					pieces[format_idx] = 'GT:' + FORMAT +':AD:AF:SBS'
					pieces[normal_idx] = './.:'+ NORMAL + ':' + str(NORMAL_AD_REF) + ',' + str(NORMAL_AD_ALT) + ':' + str(NORMAL_AF) + ':' + str(snsvs3[1])
					pieces[tumor_idx] = './.:' + TUMOR + ':' + str(TUMOR_AD_REF) + ',' + str(TUMOR_AD_ALT) + ':' + str(TUMOR_AF) + ':' + str(snsvs3[1])
					writer.write('\t'.join(pieces) + '\n')
	writer.close()
	reader.close()


def main():
# The main argument parser
	parser = ArgumentParser(description="Script for modifying Strelka Indels VCF by adding alelic depth and fraction to it")
	parser.add_argument('IN_VCF', help='input Strelka VCF containing only indels')
	parser.add_argument('OUT_VCF',help="modified output Strelka VCF ")
	parser.add_argument('TUMOR_NAME',help="tumor sample name, string that will replace TUMOR in the output Strelka VCF ")
	parser.add_argument('NORMAL_NAME',help="normal sample name, string that will replace NORMAL in the output Strelka VCF ")
	parser.add_argument("-VC", "--VariantCaller",help="Method for variant calling, will add output to INFO under ME in VCF. Strelka1: S1, Strelka2: S2", type=str, default="None")
	args = parser.parse_args()

	if(args):
		_process_vcf_file(args.IN_VCF, args.OUT_VCF, args.TUMOR_NAME, args.NORMAL_NAME, args.VariantCaller)
	else:
		parser.print_help()


if __name__ == "__main__":
	main()