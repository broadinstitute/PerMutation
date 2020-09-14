#!/usr/bin/env python
# Rauf Salamzade
# HopCountsBroad.py

import os
import sys
import re
import subprocess
import argparse
from collections import defaultdict
from Bio import SeqIO
import pysam

try:
    assert (sys.version_info[0] == 3)
except:
    sys.stderr.write("Please use Python-3 to run this program. Exiting now ...\n"); sys.exit(1)

# number of cores/threads to use for trimmomatic + bowtie2
threads = 1

def fastq_filter(inFastq, outputDir, uncompressFlag):
    # Use the fastx toolkit to clip the data and the Trimmomatic software to QC filter the reads.

    clippedFastq = outputDir + '.'.join(inFastq.split('/')[-1].split('.')[:-1]) + '.clipped.fastq.gz'
    trimmedFastq = outputDir + '.'.join(inFastq.split('/')[-1].split('.')[:-1]) + '.clipped.trimmed.fastq.gz'
    clippedQC = outputDir + "fastxClipperQC/"
    trimmedQC = outputDir + "trimmomaticQC/"
    createdir_cmd = ["mkdir", clippedQC, trimmedQC]
    subprocess.call(createdir_cmd)

    # Here I change the minimum length of a read to be maintained to 15, as k-mers of size 17-19 are 99%
    # likely to be unique to a region of a given bacterial reference according to scientific literature.
    # Additonally, I allow for reads with ambiguous 'N' characters to be retained.
    print('Clipping adapter string from reads ... ')
    fastx_cmd = ['fastx_clipper', '-o', clippedFastq, '-a', 'CCCCCCCCCCCCCCCCCCCCCCCCCCCCCC',
                 '-l15', '-n', '-Q33', '-z', '-i', inFastq]
    clipper_fastqc_cmd = ['fastqc', clippedFastq, '-o', clippedQC]
    subprocess.call(fastx_cmd)
    subprocess.call(clipper_fastqc_cmd)
    print('Finished')

    # NOTE: I got rid of a step to trim reads to 25 bp. To find the corresponding insertion site for each read
    # I employ a different algorithm.

    # Run Trimmomatic to filter reads based on Phred scores.
    print('Filtering for quality... ')
    trimmomatic_cmd = ['trimmomatic',
                       'SE', '-threads', str(threads), '-phred33', clippedFastq, trimmedFastq, 'LEADING:20',
                       'TRAILING:30',
                       'SLIDINGWINDOW:4:20', 'MINLEN:15']
    trim_fastqc_cmd = ["fastqc", trimmedFastq, "-o", trimmedQC]
    subprocess.call(trimmomatic_cmd)
    subprocess.call(trim_fastqc_cmd)
    print('Finished')
    cleanup_cmd = ["rm", clippedFastq]
    if uncompressFlag:
        cleanup_cmd += [inFastq]
    subprocess.call(cleanup_cmd)
    return (trimmedFastq)


def run_bowtie(inFastq, indexFile, outputDir, stringency):
    # Original bowtie v1 command run on Tufts servers was:
    # bowtie -q -p 4 -S -n 2 -e 70 -l 28 --maxbts 800 -k 1 -m 1 --best --phred33-quals S.aureus
    # Sample_MHcontrol2.R1.fastq.clp.trm.qc > Sample_MHcontrol2.R1.fastq.sam
    # Here I switch over to bowtie2:

    stats_file_stdout = outputDir + "bowtie2.stdout"
    stats_file_stderr = outputDir + "bowtie2.stderr"
    sam_file = outputDir + '.'.join(inFastq.split('/')[-1].split('.')[:-4]) + '.sam'
    filtered_sam_file = outputDir + '.'.join(inFastq.split('/')[-1].split('.')[:-4]) + '.filt.sam'
    bam_file = outputDir + '.'.join(inFastq.split('/')[-1].split('.')[:-4]) + '.filt.bam'
    bam_file_sorted = outputDir + '.'.join(inFastq.split('/')[-1].split('.')[:-4]) + '.filt.sorted.bam'

    # Run bowtie2 with very sensitive settings and don't change default reporting
    print('Running bowtie2 alignment of reads to reference index...')
    bowtie2_cmd = ['bowtie2', '-R', '10', '--met-file', stats_file_stdout, '-p', str(threads), '--phred33',
                   '--very-sensitive', '-x', indexFile, '-U', inFastq, '-S', sam_file, '&>', stats_file_stderr]
    subprocess.call(bowtie2_cmd)
    print('Finished alignment')

    outf = open(filtered_sam_file, 'w')
    with open(sam_file) as osf:
        for line in osf:
            if line.startswith("@"): outf.write(line)
            else:
                if stringency.lower().strip() == 'true-uni':
                    if not 'XS:i:' in line:
                        ls = line.rstrip('\n').split('\t')
                        if (ls[1] == "0" or ls[1] == "16") and (int(ls[4]) >= 3):
                            for l in ls:
                                if l.startswith("XM:i:"):
                                    mismatches = l.split('XM:i:')[1]
                                    if mismatches == "0" or mismatches == "1":
                                        outf.write(line)
                elif stringency.lower().strip() == 'true-multi':
                    ls = line.rstrip('\n').split('\t')
                    AS = None; XS = None
                    for v in line.split():
                        if v.startswith("AS:i:"):
                            AS = v.split('AS:i:')[1]
                        elif v.startswith("XS:i:"):
                            XS = v.split('XS:i:')[1]
                    if ls[1] != "4" and (ls[4] == "0" or ls[4] == "1") and (AS and AS == XS):
                        outf.write(line)
                elif stringency.lower().strip() == 'recommended':
                    ls = line.rstrip('\n').split('\t')
                    if (ls[1] == "0" or ls[1] == "16") and (int(ls[4]) >= 10):
                        for l in ls:
                            if l.startswith("XM:i:"):
                                mismatches = l.split('XM:i:')[1]
                                if mismatches == "0" or mismatches == "1":
                                    outf.write(line)
                else:
                    sys.stderr.write("The stringency specified for the mapping was not an accepted mode. Please read the help statement for the program! Exiting now ...\n")
                    raise RuntimeError
    outf.close()

    print('Filters applied!')
    sam_to_bam_cmd = ["samtools", "view", "-Sb", filtered_sam_file]
    ob = open(bam_file, 'w')
    subprocess.call(sam_to_bam_cmd, stdout=ob)
    ob.close()
    sort_bam_cmd = ["samtools", "sort", bam_file, '.'.join(bam_file_sorted.split('.')[:-1])]
    subprocess.call(sort_bam_cmd)
    index_bam_cmd = ["samtools", "index", bam_file_sorted]
    subprocess.call(index_bam_cmd)
    cleanup_cmd = ["rm", inFastq, sam_file, bam_file, filtered_sam_file]
    subprocess.call(cleanup_cmd)
    return bam_file_sorted


def hopcount_calc(bamFile, outputDir, refGenomeFile, resultName):
    # Parse FASTA Genome
    reference_fasta = {}
    ta_sites = {}
    with open(refGenomeFile) as rGF:
        for rec in SeqIO.parse(rGF, 'fasta'):
            reference_fasta[rec.id] = str(rec.seq)
            ta_sites[rec.id] = sorted([x + 1 for x in [m.start() for m in re.finditer('TA', str(rec.seq))]])

    scaffolds = reference_fasta.keys()

    print(scaffolds)
    pysamfile = pysam.AlignmentFile(bamFile, "rb")
    for i, scaffold in enumerate(scaffolds):
        plusData = defaultdict(int)
        minusData = defaultdict(int)
        for read in pysamfile.fetch(scaffold):
            direction = read.flag
            location = read.reference_start + 1
            readlength = read.reference_length

            # find insertion TA site. If insertion is on plus strand, find the closest TA site to the left of the
            # alignment reported location. Otherwise, if insertion is on the minus strand, find the closet TA site
            # to the right of the alignment reported location.
            location = int(location)
            if direction == 0:
                insertion_location = location - 2
                plusData[insertion_location] += 1
            else:
                insertion_location = location + readlength
                minusData[insertion_location] += 1

        # write results to output file!
        outf_hopcounts = open(outputDir + resultName + ".tsv", 'a+')
        outf_wigfile = open(outputDir + resultName + ".wig", 'a+')
        outf_wigfile.write('variableStep chrom=' + scaffold + '\n')
        if i == 0:
            outf_hopcounts.write(
                '\t'.join(['scaffold', 'insertion_location', 'plus_strand_counts', 'minus_strand_counts']) + '\n')
        for loc in ta_sites[scaffold]:
            outf_wigfile.write('\t'.join([str(loc), str(plusData[loc] + minusData[loc])]) + '\n')
            outf_hopcounts.write('\t'.join([scaffold, str(loc), str(plusData[loc]), str(minusData[loc])]) + '\n')
        outf_hopcounts.close()
        outf_wigfile.close()


def run_pipeline(fastqFile, indexName, refGenomeFile, outputDir, resultName, stringency, cores):
    global threads, strict_alignment
    if cores != 1:
        threads = cores

    """
    Main function that runs pipeline by calling all three functions.
    """

    fastqFile = os.path.abspath(fastqFile)
    indexName = os.path.abspath(indexName)
    refGenomeFile = os.path.abspath(refGenomeFile)
    outputDir = os.path.abspath(outputDir) + '/'


    if not os.path.isdir(outputDir):
        os.system('mkdir ' + outputDir)
    else:
        sys.stderr.write('Warning: output directory already exists.\n')

    # Uncompress fastqFile if necessary
    uncompress_flag = False
    if fastqFile.endswith('.gz'):
        fastqFile_cp = outputDir + fastqFile.split('/')[-1]
        os.system("cp " + fastqFile + " " + fastqFile_cp)
        os.system("gunzip " + fastqFile_cp)
        fastqFile = outputDir + '.'.join(fastqFile.split('/')[-1].split('.')[:-1])
        uncompress_flag = True

    # Call the fastq clip, trim, filter set
    processedFastq = fastq_filter(fastqFile, outputDir, uncompress_flag)

    # Call bowtie2 to do the mapping to the reference genomre
    bamFile = run_bowtie(processedFastq, indexName, outputDir, stringency)

    # Calculate TA insertion site loci and read counts for them
    hopcount_calc(bamFile, outputDir, refGenomeFile, resultName)

    print("Run has completed successfully!")

if __name__ == '__main__':
    # Parse arguments.

    parser = argparse.ArgumentParser(description="""
    Pipeline in Python to QC process TnSeq sequencing data, map them to reference with bowtie2, and compute \"Hop Counts\", 
    the count of insertions at each TA site.
    """)
    parser.add_argument('-f', '--fastq_file', help='location of input fastq file (has to be uncompressed).',
                        required=True)
    parser.add_argument('-i', '--index_name', help='name of the bowtie index', required=True)
    parser.add_argument('-r', '--ref_genome_file', help='fasta file for the reference genome', required=True)
    parser.add_argument('-o', '--output_dir', help='path to the output directory', required=True)
    parser.add_argument('-s', '--stringency',
                        help='specify alignment stringency, options include: true-uni, true-multi, and recommended. Default = recommended',
                        required=False, default="recommended")
    parser.add_argument('-n', '--result_name', help='prefix for final result naming', required=False, default='result')
    parser.add_argument('-c', '--cores', type=int,
                        help='number of cores to allocate per run, will be used for only trimmomatic + bowtie2 steps.',
                        required=False, default=1)
    args = parser.parse_args()

    run_pipeline(''.join(args.fastq_file), ''.join(args.index_name), ''.join(args.ref_genome_file),
                 ''.join(args.output_dir), ''.join(args.result_name), ''.join(args.stringency), args.cores)
