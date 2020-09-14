#!/usr/bin/env python
# Rauf Salamzade
# ModifiedDvalCalculator.py

import os
import sys
from sys import stderr
import argparse
from collections import defaultdict
from Bio import SeqIO
import math
import re

try: assert(sys.version_info[0] == 3)
except: stderr.write("Please use Python-3 to run this program. Exiting now ...\n"); sys.exit(1)

def get_scaffold_lengths(fasta, masked_TA_sites):

	masked_sites = defaultdict(set)
	if masked_TA_sites:
		with open(masked_TA_sites) as omts:
			for line in omts:
				line = line.rstrip('\n')
				ls = line.split('\t')
				masked_sites[ls[0]].add(int(ls[1]))

	scaf_lens = {}; scaf_tas = {}
	with open(fasta) as of:
		for rec in SeqIO.parse(of, 'fasta'):
			scaf_lens[rec.id] = len(str(rec.seq))
			scaf_tas[rec.id] = sorted([x+1 for x in [m.start() for m in re.finditer('TA', str(rec.seq).upper())] if not (x+1) in masked_sites[rec.id]])
	return scaf_lens, scaf_tas

def parse_gff(gff, upstream, trim, scaf_tas):
	trim_edge = (1.0-trim)/2.0
	pos_to_gene = {}
	genes = {}
	exon_tas = {}
	with open(gff) as og:
		for line in og:
			line = line.rstrip('\n')
			ls = line.split('\t')
			if not line.startswith("#") and ls[2] == 'gene':
				gene = ls[-1].split('ID=')[1].split(';')[0]
				scaffold = ls[0]; start = int(ls[3]); stop = int(ls[4]); strand = ls[6]
				if strand == '+': start -= upstream
				elif strand == '-': stop += upstream
				glen = abs(start-stop)+1
				# trim start loci and round up to integer
				start = int(math.floor(start + (glen*trim_edge)))
				# trim stop loci and round down to integer
				stop = int(math.ceil(stop - (glen*trim_edge)))
				# Please note, gene length is the length of the feature
				# after trimming and upstream region addition are accounted
				# for.
				if not pos_to_gene.has_key(scaffold):
					pos_to_gene[scaffold] = defaultdict(set)
					genes[scaffold] = {}
					exon_tas[scaffold] = 0
				gene_tas = [x for x in scaf_tas[scaffold] if x >= start and x <= stop]
				exon_tas[scaffold] += len(gene_tas)
				genes[scaffold][gene] = [start, stop]
				for bp in range(start, stop+1):
					pos_to_gene[scaffold][bp].add(gene)
	return genes, pos_to_gene, exon_tas

def load_counts(wig_file):
	hopcounts = {}
	scaf = None
	with open(wig_file) as ow:
		for i, line in enumerate(ow):
			line = line.rstrip('\n')
			if line.startswith("variableStep"):
				scaf = line.split('chrom=')[1]
				hopcounts[scaf] = defaultdict(float)
			else:
				pos, insert_count = line.split('\t')
				hopcounts[scaf][int(pos)] = float(insert_count)
	return hopcounts

def aggregate_counts(scaf_lens, pos_to_gene, hopcounts, exon_only):
	gene_counts = {}
	total_counts = defaultdict(float)
	for scaf in scaf_lens:
		gene_counts[scaf] = defaultdict(float)
		for bp in range(1, scaf_lens[scaf] + 1):
			if hopcounts[scaf].has_key(bp):
				if pos_to_gene[scaf].has_key(bp):
					genes = pos_to_gene[scaf][bp]
					for g in genes:
						gene_counts[scaf][g] += hopcounts[scaf][bp]
					total_counts[scaf] += hopcounts[scaf][bp]
				elif not exon_only:
					total_counts[scaf] += hopcounts[scaf][bp]
	return gene_counts, total_counts

def compute_dvals(output, genes, exon_tas, scaf_tas, gene_counts, total_counts, trim, upstream, exon_only):
	for s in genes:
		outf = open('.'.join([output, s.replace('.','-'), 'txt']), 'w')
		outf.write('gene\tscaffold\tstart\tstop\tgene_length\ttrim\tupstream\texon_only\tobserved_insertions\tdval\n')
		for g in genes[s]:
			observed_count = gene_counts[s][g]
			gstart, gend = genes[s][g]
			gene_tas = len([x for x in scaf_tas[s] if x >= gstart and x <= gend])
			expected_count = None
			if exon_only: expected_count = (float(gene_tas)/exon_tas[s])*total_counts[s]
			else: expected_count = (float(gene_tas)/len(scaf_tas[s]))*total_counts[s]
			dval = 'nan' 
			if expected_count > 0.0: dval = float(observed_count)/expected_count
			outf.write('\t'.join([str(x) for x in [g, s, gstart, gend, gene_tas, trim, upstream, exon_only, observed_count, dval]])+ '\n')
		outf.close()

def calc_dvalue_pipeline(wig_file, output, masked_TA_sites, fasta, gff, exon_only, upstream, trim):
	""" Main function which wraps everything together """
	try: assert(upstream >= 0)
	except: sys.stderr.write('Error: -u parameter must be provided positive integer or zero. Exiting now ...\n'); sys.exit(1)
	try: assert(trim > 0.0 and trim <= 1.0)
	except: sys.stderr.write('Error: -t trim must be provided positive float greater than 0.0 and less than or equal to 1.0. Exiting now ...\n'); sys.exit(1)
	try: assert(os.path.isfile(fasta) and os.path.isfile(wig_file) and os.path.isfile(gff))
	except:	sys.stderr.write('Error: Unable to locate hopcounts, gff, and/or fasta input files! Please check input and try again. Exiting now ...\n'); sys.exit(1)
	if trim != 1.0 and upstream != 0: sys.stderr.write('Warning: Using -u and -t parameters simultaneously! Please consider if this is truly what you want. Continuing ...\n')

	scaf_lens, scaf_tas = get_scaffold_lengths(fasta, masked_TA_sites)
	genes, pos_to_gene, exon_tas = parse_gff(gff, upstream, trim, scaf_tas)
	pos_counts = load_counts(wig_file)
	gene_counts, tot_counts = aggregate_counts(scaf_lens, pos_to_gene, pos_counts, exon_only)
	compute_dvals(output, genes, exon_tas, scaf_tas, gene_counts, tot_counts, trim, upstream, exon_only)

if __name__ == '__main__':
	#Parse arguments.
	parser = argparse.ArgumentParser(description="""
	Progam to compute the modified dVal for each gene, a simple metric which simultaneously normalizes for sequencing depth and TA count (instead of gene length) and quantifies fitness.
	""")

	parser.add_argument('-i', '--wig_file', help='Location of wig file with insertion information.', required=True)
	parser.add_argument('-o', '--output', help='Location and prefix of output file with Dval stats from aggregated counts.', required=True)
	parser.add_argument('-m', '--masked_TA_sites', help='File with coordinates of masked TA sites.', required=False)
	parser.add_argument('-f', '--fasta', help='Fasta file of genome assembly corresponding to GFF', required=True)
	parser.add_argument('-g', '--gff', help='GFF file from Vesper/Calhoun. Gene ID needs to be in aliases identified by the key "ID".', required=True)
	parser.add_argument('-e', '--exon_only', action='store_true', help='Calculate Dval accounting for only insertions which lie on genes.', required=False, default=False)
	parser.add_argument('-t', '--trim', type=float, help='Middle proportion of genes to consider. For instance, 0.8 specifies that the first and last 10 percent of bases in the gene will be ignored.', required=False, default=1.0) 
	parser.add_argument('-u', '--upstream', type=int, help='Include specified number of bases upstream as part of gene to test if including promotor regions maintains essentiallity. Should generally not be used in combination with trim. Default is 0.', default=0, required=False)

	args = parser.parse_args()

	calc_dvalue_pipeline(args.wig_file, args.output, args.masked_TA_sites, args.fasta, args.gff, args.exon_only, args.upstream, args.trim)