#!/usr/bin/env python
# -*- coding: utf-8 -*-
# https://github.com/shenwei356/bio_scripts
# Author     : Wei Shen
# Contact    : shenwei356@gmail.com
# LastUpdate : 2015-07-17

from __future__ import print_function, division
import argparse
import os
import sys
import gzip
from collections import defaultdict, Counter
import pickle
import numpy as np

parser = argparse.ArgumentParser(description="gff frame start coverage",
                                 epilog="https://github.com/shenwei356/bio_scripts")

parser.add_argument('genome_size_file', type=str, help='genome size file. two fields (chr and size) per line. ')
parser.add_argument('gff_file', type=str, help='gff/gtf file')
parser.add_argument('-w', '--window', type=int, default=1000, help='windows size [1000]')
parser.add_argument('-s', '--step', type=int, default=30, help='step size [30]')

args = parser.parse_args()

# read genome size file
sys.stderr.write('read genome size\n')
genomesizes = defaultdict(int)
with gzip.open(args.genome_size_file) if args.genome_size_file.endswith('.gz') else open(args.genome_size_file) as fh:
    for line in fh:
        if line.isspace() or line[0] == '#':
            continue
        data = line.rstrip().split()
        if len(data) < 2:
            sys.stderr.write('number of columns < 2! {}'.format(line))
            continue
        chr, size = data[0], data[1]
        genomesizes[chr] = int(size)

# read gff file
sys.stderr.write('read gff file\n')
coverages = defaultdict(dict)
file_cov_pickle = '{}.cov.pickle'.format(args.gff_file)
if not (os.path.exists(file_cov_pickle) and os.path.getsize(file_cov_pickle) > 0):
    with gzip.open(args.gff_file) if args.gff_file.endswith('.gz') else open(args.gff_file) as fh:
        chr = ''
        for line in fh:
            if line.isspace() or line[0] == '#':
                continue

            data = line.rstrip().split('\t')
            if len(data) != 9:
                sys.stderr.write('number of columns != 9: {}'.format(line))

            g, start, end, strand = data[0], int(data[3]), int(data[4]), data[6]
            if g != chr:
                chr = g
                coverages[chr]['+'] = np.zeros(genomesizes[chr] + 1, dtype=np.uint32)
                coverages[chr]['-'] = np.zeros(genomesizes[chr] + 1, dtype=np.uint32)
                sys.stderr.write('read chr {}\n'.format(chr))
            if strand == '+':
                coverages[chr][strand][start] += 1
                # print(chr, strand, start)
            else:
                coverages[chr][strand][end] += 1
    with open(file_cov_pickle, 'wb') as fh:
        pickle.dump(coverages, fh, pickle.HIGHEST_PROTOCOL)
else:
    with open(file_cov_pickle, 'rb') as fh:
        coverages = pickle.load(fh)


def mean_coverage(data):
    return round(sum((c for j, c in data)) / len(data), 2) if len(data) > 0 else 0

# counting
sys.stderr.write('statistics...\n')
sys.stdout.write('#{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format('chr', 'strand', 'start', 'end', 'cnt',
                                                                        'cnt_f0', 'cnt_f1', 'cnt_f2',
                                                                        'a_cov_f0', 'a_cov_f1', 'a_cov_f2'))
for chr in sorted(coverages.keys()):
    for strand in ['+', '-']:
        coverage = coverages[chr][strand]
        # for i, c in enumerate(coverage):
        #     sys.stdout.write('{}\t{}\t{}\t{}\n'.format(chr, strand, i, c))
        _end = genomesizes[chr] - args.window + 1 if genomesizes[chr] > args.window else 1
        # print(chr, strand, genomesizes[chr], _end)
        for i in np.arange(1, _end + 1, args.step, dtype=np.uint32):
            data = [(j, coverage[j]) for j in np.arange(i, i + args.window) if coverage[j] > 0]
            data_f0 = [(j, c) for j, c in data if j % 3 == 1]
            data_f1 = [(j, c) for j, c in data if j % 3 == 2]
            data_f2 = [(j, c) for j, c in data if j % 3 == 0]
            sys.stdout.write(
                '{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(chr, strand, i, i + args.window - 1, len(data),
                                                                      len(data_f0), len(data_f1), len(data_f2),
                                                                      mean_coverage(data_f0), mean_coverage(data_f1),
                                                                      mean_coverage(data_f2)))
