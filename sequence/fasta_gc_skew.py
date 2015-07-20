#!/usr/bin/env python
# https://github.com/shenwei356/bio_scripts
from __future__ import division
import sys
import argparse
from Bio import SeqIO
import numpy as np


def parse_args():
    parser = argparse.ArgumentParser(description="GC Skew",
                                     epilog="https://github.com/shenwei356/bio_scripts")

    parser.add_argument('infile', type=str, help='fasta file')
    parser.add_argument('-w', '--window', type=int, default=10000, help='window size [10000]')
    parser.add_argument('-s', '--step', type=int, default=200, help='step size [200]')
    parser.add_argument('-c', '--circular', action='store_true', help='circular genome')

    args = parser.parse_args()
    return args


def GC_Skew(seq, window=10000, step=200, circular=False):
    length, cnt = len(seq), 0
    if circular:
        end = length - step if length > step else 0
    else:
        end = length - window if length > window else 0
    locs = range(0, end + 1, step)
    GC, skew = np.zeros(len(locs)), np.zeros(len(locs))
    for i in locs:
        if i >= length - window:
            s = '{}{}'.format(seq[i:length], seq[0:window - (length - i)])
        else:
            s = seq[i:i + window]
        g, c = s.count('g') + s.count('G'), s.count('c') + s.count('C')
        GC[cnt] = (g + c) / window
        skew[cnt] = (g - c) / (g + c)
        cnt += 1
    return GC, skew


if __name__ == '__main__':
    args = parse_args()

    with open(args.infile) as fh:
        sys.stdout.write('{}\t{}\t{}\t{}\t{}\n'.format('chr', 'loc', 'gc', 'gcskew', 'accum_gcskew'))
        for seq in SeqIO.parse(fh, 'fasta'):
            sys.stderr.write('compute gcskew: {}\n'.format(seq.id))
            GC, gcskew = GC_Skew(seq.seq, window=args.window, step=args.step, circular=args.circular)
            acc = 0
            for i in range(0, len(GC)):
                gc, skew = GC[i], gcskew[i]
                acc += skew
                sys.stdout.write('{}\t{}\t{:.4f}\t{:.4f}\t{:.4f}\n'.format(seq.id, i * args.step + 1, gc, skew, acc))
