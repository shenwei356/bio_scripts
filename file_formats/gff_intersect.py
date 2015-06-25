#!/usr/bin/env python
# https://github.com/shenwei356/bio_scripts

from __future__ import print_function, division
import argparse
import sys
from collections import defaultdict, Counter
import operator
from bx.intervals.intersection import Intersecter, Interval

parser = argparse.ArgumentParser(description="gff intersect",
                                 epilog="https://github.com/shenwei356/bio_scripts")

parser.add_argument('query', type=str, help='gff file b (query)')
parser.add_argument('subject', type=str, help='gff file a (subject)')

parser.add_argument("-v", "--verbose", help='verbosely print information',
                    action="count", default=0)

args = parser.parse_args()

sys.stderr.write('building tree\n')
trees = dict()
with open(args.subject) as fh:
    genome = ''
    for line in fh:
        data = line.rstrip().split('\t')
        g, start, end, strand, product = data[0], int(data[3]), int(data[4]), data[6], data[8]
        if g != genome:
            genome = g
            trees[genome] = Intersecter()

        if strand == '-':  # complement strand
            start, end = -end, -start
        trees[genome].add_interval(Interval(start, end, value=product))

sys.stderr.write('querying\n')
with open(args.query) as fh:
    for line in fh:
        data = line.rstrip().split('\t')
        genome, start, end, strand, product = data[0], int(data[3]), int(data[4]), data[6], data[8]

        if genome not in trees:
            continue

        overlaps = trees[genome].find(start, end)
        if len(overlaps) == 0:
            continue

        overlap_data, stats = list(), Counter()
        for x in overlaps:
            s, e = x.start, x.end
            if s < 0:  # complement strand
                s, e = -x.end, -x.start

            overlap, t = 0, ''
            if s <= start:
                if e >= end:
                    #   start ======== end
                    #     s ------------- e
                    overlap = end - start + 1
                    t = 'cover'
                else:
                    #  start ======== end
                    #   s ------ e
                    overlap = e - start + 1
                    t = 'overlap.upstream' if strand == '+' else 'overlap.downstream'
            else:
                if e >= end:
                    #   start ======== end
                    #           s ------ e
                    overlap = end - s + 1
                    t = 'overlap.downstream' if strand == '+' else 'overlap.upstream'
                else:
                    #   start ======== end
                    #          s --- e
                    overlap = e - s + 1
                    t = 'embed'
            if strand == '+':
                frame = abs(s - start) % 3
            else:
                frame = abs(e - end) % 3

            stats[t] += 1
            overlap_data.append(
                [str(i) for i in [x.value, s, e, overlap, round(100 * overlap / (end - start + 1), 1), t, frame]])

        sys.stdout.write('>' + line)
        sys.stdout.write('summary: {}\n'.format(stats))
        sys.stdout.write('\t'.join(['name', 'start', 'end', 'overlap', 'overlap%', 'type', 'frame']) + '\n')
        for overlap in sorted(overlap_data, key=lambda o: (o[5], o[6], -float(o[4]))):
            sys.stdout.write('\t'.join(overlap) + '\n')
        sys.stdout.write('\n')