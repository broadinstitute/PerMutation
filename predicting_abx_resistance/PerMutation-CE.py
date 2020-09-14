import os
import sys
import argparse
from collections import defaultdict
import numpy
import math

def p_adjust_bh(p):
    """
    Using implementation from https://stackoverflow.com/questions/7450957/how-to-implement-rs-p-adjust-in-python/33532498#33532498
    Benjamini-Hochberg p-value correction for multiple hypothesis testing.
    """
    p = numpy.asfarray(p)
    by_descend = p.argsort()[::-1]
    by_orig = by_descend.argsort()
    steps = float(len(p)) / numpy.arange(len(p), 0, -1)
    q = numpy.minimum(1, numpy.minimum.accumulate(steps * p[by_descend]))
    return q[by_orig]

def read_in_wig_data(all_wigs):
    wigdata = []
    position = []
    for i, wf in enumerate(all_wigs):
        sample_counts = []
        with open(wf) as owf:
            for line in owf:
                line = line.strip()
                ls = line.split('\t')
                if line.startswith("variableStep"): continue
                if i == 0: position.append(int(ls[0]))
                sample_counts.append(float(ls[1]))
        wigdata.append(sample_counts)

    return([position, wigdata])

def resampling(control_data, treatment_data, S=100000):
    count_greater = 1 # to prevent p-value of 0

    n_control = len(list(control_data.values())[0])
    control_mean = numpy.average([sum(x[1]) for x in control_data.items()])
    treatment_mean = numpy.average([sum(x[1]) for x in treatment_data.items()])

    log2FC = (control_mean + 0.1)/(treatment_mean+0.1)

    stat_obs = 0
    test_data = []
    for ind in control_data:
        test_data.append(control_data[ind] + treatment_data[ind])
        control_sum = sum(control_data[ind])
        treat_sum = sum(treatment_data[ind])
        stat_obs += control_sum - treat_sum

    for j, s in enumerate(range(S)):
        stat = 0
        for pos_vals in test_data:
            perm_vals = numpy.random.permutation(pos_vals)
            stat += sum(perm_vals[:n_control])-sum(perm_vals[n_control:])
        if stat >= stat_obs: count_greater += 1
    pval_greater = float(count_greater) / float(S)

    return (stat_obs, log2FC, pval_greater)

def create_parser():
    # Parse arguments.

    parser = argparse.ArgumentParser(description=""" This program runs a position aware variant of the permutation test developed by DeJesus et al 2015.""")

    parser.add_argument('-c', '--control_wigs', help='control wig files separated by comma.', required=True)
    parser.add_argument('-e', '--experiment_wigs', help='experimental wig files separated by comma.', required=True)
    parser.add_argument('-b', '--bed_file', help='BED file with four columns: scaffold id, start position, end position, gene id', required=True)
    parser.add_argument('-s', '--scaffold', help='scaffold id', required=True)
    parser.add_argument('-o', '--output', help='output file', required=True)

    args = parser.parse_args()
    return args

# parse arguments
myargs = create_parser()
control_wigs = myargs.control_wigs.split(',')
experiment_wigs = myargs.experiment_wigs.split(',')
bed_file = myargs.bed_file
output_file = myargs.output
scaffold = myargs.scaffold

replicate_labels = ['control']*len(control_wigs) + ['experiment']*len(experiment_wigs)
all_wigs = control_wigs + experiment_wigs
position, wigdata = read_in_wig_data(all_wigs)

pos_to_index = {}
for i, p in enumerate(position):
    nz_flag = False
    for r in wigdata:
        if r[i] > 0.0:
            nz_flag = True
    if nz_flag:
        pos_to_index[p] = i

all_info = []
pvalues = []
with open(bed_file) as obf:
    for line in obf:
            line = line.rstrip('\n')
            feature_scaffold, start, stop, description = line.split('\t')
            relevant_indices = []
            start = min(int(start), int(stop))
            stop = max(int(start), int(stop))
            glen = abs(start - stop) + 1
            start = int(math.floor(start + (glen * 0.10)))
            stop = int(math.ceil(stop - (glen * 0.10)))
            for p in range(start, stop+1):
                if p in pos_to_index:
                    relevant_indices.append(pos_to_index[p])
            if not relevant_indices: continue
            control_samples = defaultdict(list)
            experiment_samples = defaultdict(list)
            for ind in relevant_indices:
                for r in range(0,len(control_wigs)):
                    control_samples[ind].append(wigdata[r][ind])
                for r in range(len(control_wigs), (len(control_wigs)+len(experiment_wigs))):
                    experiment_samples[ind].append(wigdata[r][ind])

            stat, log2FC, pval = resampling(control_samples, experiment_samples)
            all_info.append([feature_scaffold, start, stop, description, log2FC, stat, pval])
            pvalues.append(pval)

outf = open(output_file, 'w')
qvalues = p_adjust_bh(pvalues)

for i, q in enumerate(qvalues):
    outf.write('\t'.join([ str(x) for x in (all_info[i] + [q])])+ '\n')
outf.close()
