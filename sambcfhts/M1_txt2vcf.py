#/usr/local/bin/python
import sys
from argparse import ArgumentParser


def check_allele(ref, alt, check):
	if check == ref:
		al = "0"
	elif check == alt:
		al = "1"
	else:
		al = "."
	return al

info_field_dict={"t_lod_fstar":"TLOD","n_lod_fstar":"NLOD"}
info_field_dict2={"t_lod_fstar":"TLOD_M1","n_lod_fstar":"NLOD_M1"}

def main(args):
	parser = ArgumentParser(description='Convert Mutect1 txt output to a vcf')
	parser.add_argument('in_call_stats',help='input txt file')
	parser.add_argument('outFile',help='output vcf file')
	parser.add_argument('case_name',help="insert tumour sample name")
	parser.add_argument('ctrl_name',help="insert normal control name")
	parser.add_argument("-VC", "--VariantCaller", help="VariantCallerUsed", type=str, default="None")
	parser.add_argument("-RM", "--runMode", help="'Paired' or 'TumOnly' ?", type=str, default="Paired")
 	
	args = parser.parse_args()
	if(args):
		print('Picked up   in_call_stats : '+str(args.in_call_stats))
		print('Picked up   outFile : '+str(args.outFile))
		print('Picked up   case_name : '+str(args.case_name))
		print('Picked up   ctrl_name : '+str(args.ctrl_name))
		print('Picked up   methodology : '+str(args.VariantCaller))
		print('Picked up   runMode : '+str(args.runMode))		
		all_lines(args.in_call_stats, args.outFile, args.case_name, args.ctrl_name, args.VariantCaller, args.runMode)
	else:
		parser.print_help()

def all_lines(in_call_stats, outFile, case_name, ctrl_name, VC, runMode):
	with open(in_call_stats) as cs, open(outFile, "w") as headerless_vcf:
			headerless_vcf.write("""##fileformat=VCFv4.2
##FILTER=<ID=alt_allele_in_normal,Description="Evidence seen in the normal sample">
##FILTER=<ID=clustered_read_position,Description="Clustered events observed in the tumor">
##FILTER=<ID=germline_risk,Description="Evidence indicates this site is germline, not somatic">
##FILTER=<ID=homologous_mapping_event,Description="More than three events were observed in the tumor">
##FILTER=<ID=alt_allele_in_normal,Description="">
##FILTER=<ID=seen_in_panel_of_normals,Description="Seen in at least 2 samples in the panel of normals">
##FILTER=<ID=str_contraction,Description="Site filtered due to contraction of short tandem repeat region">
##FILTER=<ID=fstar_tumor_lod,Description="Tumor does not meet likelihood threshold">
##FILTER=<ID=normal_lod,Description="">
##FILTER=<ID=germline_risk,Description="">
##FILTER=<ID=poor_mapping_region_alternate_allele_mapq,Description="">
##FILTER=<ID=possible_contamination,Description="">
##FILTER=<ID=strand_artifact,Description="">
##FILTER=<ID=normal_lod,Description="">
##FILTER=<ID=poor_mapping_region_mapq0,Description="">
##FILTER=<ID=nearby_gap_events,Description="">
##FILTER=<ID=triallelic_site,Description="Site filtered because more than two alt alleles pass tumor LOD">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">
##FORMAT=<ID=AF,Number=A,Type=Float,Description="Allele fraction of the event in the tumor">
##FORMAT=<ID=ALT_F1R2,Number=1,Type=Integer,Description="Count of reads in F1R2 pair orientation supporting the alternate allele">
##FORMAT=<ID=ALT_F2R1,Number=1,Type=Integer,Description="Count of reads in F2R1 pair orientation supporting the alternate allele">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth (reads with MQ=255 or with bad mates are filtered)">
##FORMAT=<ID=FOXOG,Number=1,Type=Float,Description="Fraction of alt reads indicating OxoG error">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=PGT,Number=1,Type=String,Description="Physical phasing haplotype information, describing how the alternate alleles are phased in relation to one another">
##FORMAT=<ID=PID,Number=1,Type=String,Description="Physical phasing ID information, where each unique ID within a given sample (but not across samples) connects records within a phasing group">
##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification">
##FORMAT=<ID=QSS,Number=1,Type=Integer,Description="Sum of base quality scores for each allele">
##INFO=<ID=DB,Number=0,Type=Flag,Description="dbSNP Membership">
##INFO=<ID=COVERED,Number=0,Type=Flag,Description="">
##INFO=<ID=ECNT,Number=1,Type=Integer,Description="Number of events in this haplotype">
##INFO=<ID=HCNT,Number=1,Type=Integer,Description="Number of haplotypes that support this variant">
##INFO=<ID=MAX_ED,Number=1,Type=Integer,Description="Maximum distance between events in this active region">
##INFO=<ID=MIN_ED,Number=1,Type=Integer,Description="Minimum distance between events in this active region">
##INFO=<ID=NLOD_M1,Number=1,Type=Float,Description="Normal LOD score">
##INFO=<ID=PON,Number=1,Type=String,Description="Count from Panel of Normals">
##INFO=<ID=RPA,Number=.,Type=Integer,Description="Number of times tandem repeat unit is repeated, for each allele (including reference)">
##INFO=<ID=RU,Number=1,Type=String,Description="Tandem repeat unit (bases)">
##INFO=<ID=STR,Number=0,Type=Flag,Description="Variant is a short tandem repeat">
##INFO=<ID=TLOD_M1,Number=1,Type=Float,Description="Tumor LOD score">
##Mutect Version=1.1.6										
##contig=<ID=1,length=249250621,assembly=GRCh37>										
##contig=<ID=2,length=243199373,assembly=GRCh37>										
##contig=<ID=3,length=198022430,assembly=GRCh37>										
##contig=<ID=4,length=191154276,assembly=GRCh37>										
##contig=<ID=5,length=180915260,assembly=GRCh37>										
##contig=<ID=6,length=171115067,assembly=GRCh37>										
##contig=<ID=7,length=159138663,assembly=GRCh37>										
##contig=<ID=8,length=146364022,assembly=GRCh37>										
##contig=<ID=9,length=141213431,assembly=GRCh37>										
##contig=<ID=10,length=135534747,assembly=GRCh37>										
##contig=<ID=11,length=135006516,assembly=GRCh37>										
##contig=<ID=12,length=133851895,assembly=GRCh37>										
##contig=<ID=13,length=115169878,assembly=GRCh37>										
##contig=<ID=14,length=107349540,assembly=GRCh37>										
##contig=<ID=15,length=102531392,assembly=GRCh37>										
##contig=<ID=16,length=90354753,assembly=GRCh37>										
##contig=<ID=17,length=81195210,assembly=GRCh37>										
##contig=<ID=18,length=78077248,assembly=GRCh37>										
##contig=<ID=19,length=59128983,assembly=GRCh37>										
##contig=<ID=20,length=63025520,assembly=GRCh37>										
##contig=<ID=21,length=48129895,assembly=GRCh37>										
##contig=<ID=22,length=51304566,assembly=GRCh37>										
##contig=<ID=X,length=155270560,assembly=GRCh37>										
##contig=<ID=Y,length=59373566,assembly=GRCh37>										
##contig=<ID=MT,length=16569,assembly=GRCh37>										
##contig=<ID=GL000207.1,length=4262,assembly=GRCh37>										
##contig=<ID=GL000226.1,length=15008,assembly=GRCh37>										
##contig=<ID=GL000229.1,length=19913,assembly=GRCh37>										
##contig=<ID=GL000231.1,length=27386,assembly=GRCh37>										
##contig=<ID=GL000210.1,length=27682,assembly=GRCh37>										
##contig=<ID=GL000239.1,length=33824,assembly=GRCh37>										
##contig=<ID=GL000235.1,length=34474,assembly=GRCh37>										
##contig=<ID=GL000201.1,length=36148,assembly=GRCh37>										
##contig=<ID=GL000247.1,length=36422,assembly=GRCh37>										
##contig=<ID=GL000245.1,length=36651,assembly=GRCh37>										
##contig=<ID=GL000197.1,length=37175,assembly=GRCh37>										
##contig=<ID=GL000203.1,length=37498,assembly=GRCh37>										
##contig=<ID=GL000246.1,length=38154,assembly=GRCh37>										
##contig=<ID=GL000249.1,length=38502,assembly=GRCh37>										
##contig=<ID=GL000196.1,length=38914,assembly=GRCh37>										
##contig=<ID=GL000248.1,length=39786,assembly=GRCh37>										
##contig=<ID=GL000244.1,length=39929,assembly=GRCh37>										
##contig=<ID=GL000238.1,length=39939,assembly=GRCh37>										
##contig=<ID=GL000202.1,length=40103,assembly=GRCh37>										
##contig=<ID=GL000234.1,length=40531,assembly=GRCh37>										
##contig=<ID=GL000232.1,length=40652,assembly=GRCh37>										
##contig=<ID=GL000206.1,length=41001,assembly=GRCh37>										
##contig=<ID=GL000240.1,length=41933,assembly=GRCh37>										
##contig=<ID=GL000236.1,length=41934,assembly=GRCh37>										
##contig=<ID=GL000241.1,length=42152,assembly=GRCh37>										
##contig=<ID=GL000243.1,length=43341,assembly=GRCh37>										
##contig=<ID=GL000242.1,length=43523,assembly=GRCh37>										
##contig=<ID=GL000230.1,length=43691,assembly=GRCh37>										
##contig=<ID=GL000237.1,length=45867,assembly=GRCh37>										
##contig=<ID=GL000233.1,length=45941,assembly=GRCh37>										
##contig=<ID=GL000204.1,length=81310,assembly=GRCh37>										
##contig=<ID=GL000198.1,length=90085,assembly=GRCh37>										
##contig=<ID=GL000208.1,length=92689,assembly=GRCh37>										
##contig=<ID=GL000191.1,length=106433,assembly=GRCh37>										
##contig=<ID=GL000227.1,length=128374,assembly=GRCh37>										
##contig=<ID=GL000228.1,length=129120,assembly=GRCh37>										
##contig=<ID=GL000214.1,length=137718,assembly=GRCh37>										
##contig=<ID=GL000221.1,length=155397,assembly=GRCh37>										
##contig=<ID=GL000209.1,length=159169,assembly=GRCh37>										
##contig=<ID=GL000218.1,length=161147,assembly=GRCh37>										
##contig=<ID=GL000220.1,length=161802,assembly=GRCh37>										
##contig=<ID=GL000213.1,length=164239,assembly=GRCh37>										
##contig=<ID=GL000211.1,length=166566,assembly=GRCh37>										
##contig=<ID=GL000199.1,length=169874,assembly=GRCh37>										
##contig=<ID=GL000217.1,length=172149,assembly=GRCh37>										
##contig=<ID=GL000216.1,length=172294,assembly=GRCh37>										
##contig=<ID=GL000215.1,length=172545,assembly=GRCh37>										
##contig=<ID=GL000205.1,length=174588,assembly=GRCh37>										
##contig=<ID=GL000219.1,length=179198,assembly=GRCh37>										
##contig=<ID=GL000224.1,length=179693,assembly=GRCh37>										
##contig=<ID=GL000223.1,length=180455,assembly=GRCh37>										
##contig=<ID=GL000195.1,length=182896,assembly=GRCh37>										
##contig=<ID=GL000212.1,length=186858,assembly=GRCh37>										
##contig=<ID=GL000222.1,length=186861,assembly=GRCh37>										
##contig=<ID=GL000200.1,length=187035,assembly=GRCh37>										
##contig=<ID=GL000193.1,length=189789,assembly=GRCh37>										
##contig=<ID=GL000194.1,length=191469,assembly=GRCh37>										
##contig=<ID=GL000225.1,length=211173,assembly=GRCh37>										
##contig=<ID=GL000192.1,length=547496,assembly=GRCh37>										
##contig=<ID=NC_007605,length=171823,assembly=NC_007605.1>										
""")
			header = cs.readline()
			while header[0] == "#" or not header.strip():
				header = cs.readline()
			header = header.strip("\n").split("\t")

			h = {x[1].lower(): x[0] for x in enumerate(header)}
			h["."] = -1
			used_vals = set(["dbsnp_site", "covered", "tumor_name", "normal_name"] + ["normal_f", "n_ref_count", "n_alt_count", "n_ref_sum", "n_alt_sum"] + ["tumor_f", "t_ref_count", "t_alt_count", "t_ref_sum", "t_alt_sum"] + ["contig", "position", ".", "ref_allele", "alt_allele", "."] + ["failure_reasons", "judgement"] + ["normal_best_gt"])
			for info_field in h:
				if info_field not in used_vals and info_field not in info_field_dict:
					headerless_vcf.write('##INFO=<ID='+info_field+',Number=1,Type=String,Description="internal mutect statistic">\n')
			if runMode=="Paired":		
				headerless_vcf.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + ctrl_name + "\t" + case_name + "\n")
			else:
				headerless_vcf.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + case_name + "\n")

			for line in cs:
				vals = line.strip("\n").split("\t")
				if len(vals) == 1:
					continue  # skip blank lines
				vals.append(".")
				vcf_line = list(map(lambda x: vals[h[x]], ["contig", "position", ".", "ref_allele", "alt_allele", "."]))
				filter_val = "PASS" if vals[h["judgement"]] == "KEEP" else vals[h["failure_reasons"]].replace(",", ";")
				vcf_line.append(filter_val)
				na1, na2 = vals[h["normal_best_gt"]]
				ngt = "/".join(list(map(lambda x: check_allele(vcf_line[3], vcf_line[4], x), [na1, na2])))
				nad = ",".join(list(map(lambda x: vals[h[x]], ["n_ref_count", "n_alt_count"])))
				nqss = ",".join(list(map(lambda x: vals[h[x]], ["n_ref_sum", "n_alt_sum"])))

				normal_format = ":".join([ngt,nad,nqss] + list(map(lambda x: vals[h[x]], ["normal_f"])))

				# tgt = "1/1" if ngt=="1/1" else "./."
				tgt = "./."
				tad = ",".join(list(map(lambda x: vals[h[x]], ["t_ref_count", "t_alt_count"])))
				tqss = ",".join(list(map(lambda x: vals[h[x]], ["t_ref_sum", "t_alt_sum"])))

				tumor_format = ":".join([tgt,tad,tqss] + list(map(lambda x: vals[h[x]], ["tumor_f"])))

				info = []
				if vals[h["covered"]] == "COVERED":
					info.append("COVERED")

				if vals[h["dbsnp_site"]] == "DBSNP":
					info.append("DB")

				format_names = ["GT","AD","QSS","AF"]
				format_names = ":".join(format_names)
				for key, value in h.items():
					if key not in used_vals:
						if key not in info_field_dict:
							info.append(key + "=" + vals[value])
						else:
							info.append(info_field_dict2[key] + "=" + vals[value])
				if VC!="None":
					info.append("ME=" + VC)
				info = ";".join(info)
				vcf_line.append(info)
				vcf_line.append(format_names)
				if runMode=="Paired":
					vcf_line.append(normal_format)
				vcf_line.append(tumor_format)
				headerless_vcf.write("\t".join(vcf_line) + "\n")

#===============================================================================
# Run Main
#===============================================================================
if __name__ == '__main__':
	main(sys.argv[1:])
	