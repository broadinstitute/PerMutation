#!/usr/bin/env python
# Rauf Salamzade
# DvalCalculator.py
import os
import sys
import subprocess
from sys import stderr
import argparse
from collections import defaultdict
from Bio import SeqIO
import math
import operator

try: assert(sys.version_info[0] == 3)
except: stderr.write("Please use Python-3 to run this program. Exiting now ...\n"); sys.exit(1)

def correction_pipeline(wig_file, chrom_scaffold, out_file, zero_deflate, window_size, rscript_path):
	tmp_dir = os.path.abspath(out_file + '_tmp_' + str(os.getpid())) + '/'
	try: assert(not os.path.exists(tmp_dir))
	except: sys.stderr.write('Error: somehow the tmp directory path is already in existence. Please try running again. Exiting now ...\n'); sys.exit(1)
	else: os.system('mkdir ' + tmp_dir)
	
	scaffold = None
	outf = None
	corrected_results = defaultdict(list)
	for i, line in enumerate(open(wig_file)):
		line = line.rstrip('\n')
		if line.startswith("variableStep"):
			if i != 0:
				if scaffold == chrom_scaffold:
					outf.close()
					rscript_cmd = ['Rscript', rscript_path, tmp_dir + scaffold + '.txt', str(window_size), tmp_dir + scaffold + '.adjusted.txt']
					subprocess.call(rscript_cmd)
					for j, l in enumerate(open(tmp_dir + scaffold + '.adjusted.txt')):
						if j > 0:
							l = l.rstrip('\n')
							corrected_results[scaffold].append([float(x) for x in l.split('\t')])
			scaffold = line.split('chrom=')[1]
			if scaffold == chrom_scaffold:
				outf = open(tmp_dir + scaffold + '.txt', 'w')
				outf.write('loci\tcounts\n')
		elif scaffold == chrom_scaffold:
			if zero_deflate and float(line.split('\t')[1]) != 0.0: outf.write(line + '\n')
			elif not zero_deflate: outf.write(line + '\n')
	if i > 0:
		if scaffold == chrom_scaffold:
			outf.close()
			rscript_cmd = ['Rscript', rscript_path, tmp_dir + scaffold + '.txt', str(window_size), tmp_dir + scaffold+  '.adjusted.txt']
			subprocess.call(rscript_cmd)
			for j, l in enumerate(open(tmp_dir + scaffold + '.adjusted.txt')):
				if j > 0:
					l = l.strip('\n')
					corrected_results[scaffold].append([float(x) for x in l.split('\t')])
			
	of = open(out_file, 'w')
	for s in corrected_results:
		of.write('variableStep chrom=' + s + '\n')
		for loc in sorted(corrected_results[s], key=operator.itemgetter(0)):
			of.write(str(int(loc[0])) + '\t' + str(loc[1]) + '\n')
	of.close()

	os.system('rm -rf ' + tmp_dir)

if __name__ == '__main__':
	#Parse arguments.
	parser = argparse.ArgumentParser(description="""
	This is a wrapper to a TnseqDiff R function which applies LOESS correction for replication bias. Replication bias
	correction is only performed for the chromosome, which is assumed to have a complete representation in the assembly. 
	""")

	parser.add_argument('-i', '--wig_file', type=str, help='Location of wig file with insertion information.', required=True)
	parser.add_argument('-c', '--chrom_scaffold', type=str, help='Name of the chromosomal scaffold in the assembly/wig-file.', required=True)
	parser.add_argument('-o', '--out_file', type=str, help='Location of resulting wig file where counts have been corrected for replication bias.', required=True)
	parser.add_argument('-z', '--zero_deflate', action='store_true', help='Whether to exclude zero-insertion TA sites when correcting for replication bias. Default is false.', required=False, default=False)
	parser.add_argument('-w', '--window_size', type=int, help='The window size to use when correcting for replication bias. Default is 10,000 .', required=False, default=10000)
	parser.add_argument('-p', '--rscript_path', type=str, help='Path to Rscript ReplicationCorrector.R for replication bias correction.', required=True)
	args = parser.parse_args()

	correction_pipeline(args.wig_file, args.chrom_scaffold, args.out_file, args.zero_deflate, args.window_size, args.rscript_path)
