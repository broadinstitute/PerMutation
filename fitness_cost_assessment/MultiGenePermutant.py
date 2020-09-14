import os
import sys
import numpy as np
import argparse
import random
from scipy import stats
import math

def permutant(wig, gff, output, start, stop, window, genome_length, permutation):

    genic_regions = set([])

    with open(gff) as ogf:
        for line in ogf:
            line = line.rstrip()
            ls = line.split('\t')
            if ls[2] != 'gene': continue
            for p in range(int(ls[3]), int(ls[4])+1): genic_regions.add(p)

    window_start = [start - window, None]
    window_stop = [stop + window, None]
    if window_start[0] < 1:
        alt_start = genome_length - (window - start)
        window_start = [1, alt_start]
        window_stop[1] = genome_length
    elif window_stop[0] > genome_length:
        alt_stop = window_stop[0] - genome_length
        window_stop = [genome_length, alt_stop]
        window_start[1] = 1

    window_tot_ta = 0
    window_sat_ta = 0
    region_ta = 0
    region_sum = 0
    window_nz_counts = []
    continuous_zeros = 0
    longest_continuous_zeros = 0
    with open(wig) as ow:
        for line in ow:
            line = line.rstrip('\n')
            if line.startswith('variable'): continue
            ls = line.split('\t')
            pos = int(ls[0])
            count = float(ls[1])
            if ((pos >= window_start[0] and pos <= window_stop[0]) or (window_start[1] and pos >= window_start[1] and pos <= window_stop[1])) and pos in genic_regions:
                if pos >= start and pos <= stop:
                    region_ta += 1
                    region_sum += count
                    if count > 0.0:
                        if continuous_zeros > longest_continuous_zeros: longest_continuous_zeros = continuous_zeros
                        continuous_zeros = 0
                    else:
                        continuous_zeros += 1
                window_tot_ta += 1
                if count > 0.0:
                    window_nz_counts.append(count)
                    window_sat_ta += 1

    if continuous_zeros > longest_continuous_zeros: longest_continuous_zeros = continuous_zeros
    p = window_sat_ta/float(window_tot_ta)
    empirical_pvalue = 0

    for i in range(0, permutation):
        simulated_nz_sites = np.random.binomial(region_ta, p, 1)[0]
        simulated_sum = 0
        random.shuffle(window_nz_counts)
        if simulated_nz_sites > 0: simulated_sum = sum(window_nz_counts[:simulated_nz_sites])
        if simulated_sum <= region_sum: empirical_pvalue += 1

    outf = open(output, 'w')
    outf.write(str(empirical_pvalue/permutation) + '\t' + str(p) + '\n')
    outf.close()

if __name__ == '__main__':
    # Pull out the arguments.
    parser = argparse.ArgumentParser(description="""
	Progam to compute the Dval for each gene, a simple metric which simultaneously normalizes for sequencing depth and gene length and quantifies fitness.
	""")

    parser.add_argument('-i', '--wig_file', help='Location of wig file with insertion information.', required=True)
    parser.add_argument('-a', '--gff', help='gff annotation', required=True)
    parser.add_argument('-o', '--output', help='Location and prefix of output file with Dval stats from aggregated counts.', required=True)
    parser.add_argument('-s', '--start', type=int, help='start position', required=True)
    parser.add_argument('-t', '--stop', type=int, help='stop position', required=True)
    parser.add_argument('-w', '--window', type=int, help='The window size.', required=True)
    parser.add_argument('-g', '--genome_length', type=int, help='Length of genome', required=True)
    parser.add_argument('-p', '--permutation', type=int, help='permutations', required=False, default=100000)

    args = parser.parse_args()
    permutant(args.wig_file, args.gff, args.output, args.start, args.stop, args.window, args.genome_length, args.permutation)


