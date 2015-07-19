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
    parser.add_argument('-w', '--window', type=int, default=10000, help='window size')
    parser.add_argument('-s', '--step', type=int, default=200, help='step size')

    args = parser.parse_args()
    return args


def GC_skew(seq, window=10000, step=200):
    length, cnt = len(seq), 0
    end = length - window if length > window else 0
    locs = range(0, end + 1, step)
    skew = np.zeros(len(locs))
    for i in locs:
        s = seq[i:i + window]
        g, c = s.count('g') + s.count('G'), s.count('c') + s.count('C')
        skew[cnt] = (g - c) / (g + c)
        cnt += 1
    return skew


if __name__ == '__main__':
    args = parse_args()

    with open(args.infile) as fh:
        sys.stdout.write('{}\t{}\t{}\t{}\n'.format('chr', 'loc', 'gcskew', 'accum_gcskew'))
        for seq in SeqIO.parse(fh, 'fasta'):
            sys.stderr.write('compute gcskew: {}\n'.format(seq.id))
            gcskew = GC_skew(seq.seq, window=args.window, step=args.step)
            acc = 0
            for i, skew in enumerate(gcskew):
                acc += skew
                sys.stdout.write('{}\t{}\t{:.4f}\t{:.4f}\n'.format(seq.id, i*args.step+1, skew, acc))
