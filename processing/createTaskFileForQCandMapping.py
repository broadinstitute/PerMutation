import os
import sys
import argparse

def create_parser():
    # Parse arguments.

    parser = argparse.ArgumentParser(description="""
		This quick script will generate a "task" file where each line corresponds to a HopCountsBroad.py command for
		processing and mapping a single fastq/sample file.""")

    parser.add_argument('-f', '--fastq_dir', help='location of directory with input fastq files.', required=True)
    parser.add_argument('-i', '--index_name', help='name of the bowtie index', required=True)
    parser.add_argument('-r', '--ref_genome_file', help='fasta file for the reference genome', required=True)
    parser.add_argument('-o', '--output_dir', help='path to the general output directory', required=True)
    parser.add_argument('-p', '--program', help='path to Python program HopCountsBroad.py.', required=True)
    parser.add_argument('-c', '--cores', type=int,
                        help='number of cores to allocate per run, will be used for only trimmomatic + bowtie2 steps.',
                        required=False, default=1)
    parser.add_argument('-m', '--mapping', help="Provide mapping of sequencing data to sample names. First column is sample name, second column is sequencing file [up to .fastq or .fastq.gz]",
                        required=False, default=None)
    parser.add_argument('-s', '--stringency', help='specify which mapping stringency to use.', required=False, default="recommended")
    args = parser.parse_args()

    return ''.join(args.fastq_dir), ''.join(args.index_name), ''.join(args.ref_genome_file), ''.join(args.output_dir), args.cores, ''.join(args.mapping), ''.join(args.stringency), ''.join(args.program)

fastq_dir, index_name, ref_genome_file, output_dir, cores, mapping_file, stringency, program = create_parser()

fastq_dir = os.path.abspath(fastq_dir) + '/'
output_dir = os.path.abspath(output_dir) + '/'
ref_genome_file = os.path.abspath(ref_genome_file)
index_name = os.path.abspath(index_name)
program = os.path.abspath(program)
mapping_file = os.path.abspath(mapping_file)
stringency = stringency

mapping = {}
if os.path.isfile(mapping_file):
    with open(mapping_file) as omf:
        for line in omf:
            line = line.rstrip('\n')
            ls = line.split('\t')
            mapping[ls[1].split('.fastq')[0]] = ls[0]

if not os.path.isdir(output_dir):
    os.system('mkdir ' + output_dir)

for f in os.listdir(fastq_dir):
    if f.endswith(".fastq") or f.endswith(".fastq.gz"):
        fastqFile = fastq_dir + f
        outDir = f.split('.fastq')[0]
        if not f.startswith("Undetermined") and outDir in mapping:
            if mapping:
                outDir = mapping[outDir]
            print(' '.join([program, '-f', fastqFile, '-i', index_name, '-r', ref_genome_file, '-o',
                            output_dir + outDir, '-n', outDir, '-c', str(cores), '-s', stringency]))
