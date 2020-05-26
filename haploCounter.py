#!/usr/bin/env python

"""
Author: Geert van Geest
"""

import sys
import argparse
import pysam
import numpy as np

def vcf_to_list(vcffile, regiontup):
	'''
	for a region get the positions of the SNPs
	'''
	bcf_in = pysam.VariantFile(vcffile)  # auto-detect input format
	snpl = []
	for rec in bcf_in.fetch(regiontup[0], regiontup[1], regiontup[2]):
		snpl.append(rec.pos)
	return snpl	


def frg_from_pileup(snpl, samfile, regiontup):
	'''
	from snp positions in a region get the alleles for each read from a samfile
	'''
	samfile = pysam.AlignmentFile(samfile, "rb")
	nucl_list = []
	nucl_list_agg = []
	read_list_agg = []
	for loc in snpl:
		ploc = loc - 1
		nl = []
		rl = []
		for pileupcolumn in samfile.pileup(regiontup[0], ploc -2, ploc +2):
			pileupcolumn.set_min_base_quality(0)
			if pileupcolumn.pos == ploc:
				for pileupread in pileupcolumn.pileups:
					rl.append(pileupread.alignment.query_name)
					if pileupread.is_del or pileupread.is_refskip:
						nl.append("-")
					else:
						nl.append(pileupread.alignment.query_sequence[pileupread.query_position])	
		nucl_list_agg.append(nl)
		read_list_agg.append(rl)
		
	nucl_arr = np.array(nucl_list_agg)

	nucl_arrt = np.transpose(nucl_arr)
	
	nj_list = []
	for i in range(len(nucl_arrt)):
		nj_list.append("".join(nucl_arrt[i,]))
	
	if(len(nj_list)>0):
		return([read_list_agg[0], nj_list])
	else:
		return(None)

def count_haplo(haplolist):
	'''
	count the number of individual haplotypes from a list of haplotypes
	'''
	if haplolist is None:
		return(None)
	hl = haplolist[1]
	hlu = set(hl)
	hcount_dict = {}
	for h in hlu:
		hcount_dict[h] = hl.count(h)
	return(hcount_dict)


def write_counts(count_dict, sample, regtup, con, mincount):
	'''
	write output to tab-delimited file
	'''
	if count_dict is None:
		return(None)
	#con = open(connection, 'a')
	for hap in count_dict.keys():
		if count_dict[hap] >= mincount:
			out = [sample, regtup[0], str(regtup[1]), str(regtup[2]), hap, str(count_dict[hap])]	
			con.write("	".join(out) + "\n")
	#con.close()
	
def run_haplocounter(bed, vcf, sam, sample, outfile, mincount):
	'''
	loop over bed file and run separate functions
	'''
	bedcon = open(bed, "r")
	outcon = open(outfile, "w")
	for line in bedcon:
		reg = line.split("\t")
		regtup = (reg[0], int(reg[1]), int(reg[2]))
		snpl = vcf_to_list(vcf, regtup)
		haplolist = frg_from_pileup(snpl, sam, regtup)		
		ch = count_haplo(haplolist)
		write_counts(ch, sample, regtup, outcon, mincount)
	bedcon.close()
	outcon.close()


if __name__ == "__main__":
	description_text = 'Processing of BAM files into haplotype counts.'
	epilog_text = 'Copyright GvG'
	# parser = argparse.ArgumentParser(description=description_text, epilog=epilog_text)
	# parser.add_argument('-i', type=str, required=True, dest='input.bam', help='Input BAM file. Reads should cover entire specified region in -r')
	# parser.add_argument('-s', type=str, required=True, dest='sampleID', help='Sample ID')
	# parser.add_argument('-vcf', type=str, required=True, dest='SNPS.vcf.gz', help='Input (tabix-indexed) VCF/BCF file')
	# parser.add_argument('-r', type=str, required=True, dest='CHROM:start-end', help='Region chr:start-end')
	# parser.add_argument('-c', type=int, default = 0, help='Minimum number of reads supporting haplotype')
	# parser.add_argument('-o', type=str, default='-', help='Output count file')
	# parser.add_argument('-wmode', type=str, default='w', help="Output write mode. E.g. 'a' for append or 'w' for write ")

	parser = argparse.ArgumentParser(description=description_text, epilog=epilog_text)
	parser.add_argument('-i', type=str, required=True, help='Input BAM file. Reads should cover entire specified region in -r')
	parser.add_argument('-s', type=str, required=True, help='Sample ID')
	parser.add_argument('-vcf', type=str, required=True, help='Input (tabix-indexed) VCF/BCF file')
	parser.add_argument('-bed', type=str, required=True, help='Input bed file with target regions')
	parser.add_argument('-c', type=int, default = 0, help='Minimum number of reads supporting haplotype')
	parser.add_argument('-o', type=str, default='-', help='Output count file')
	# parser.add_argument('-wmode', type=str, default='w', help="Output write mode. E.g. 'a' for append or 'w' for write ")


	args = parser.parse_args()
	
	run_haplocounter(args.bed, args.vcf, args.i, args.s, args.o, args.c)
	
	# rspl = args.r.split(":")
	# chr = rspl[0]
	# reg = rspl[1].split("-")
	# start = int(reg[0])
	# end = int(reg[1])
	
