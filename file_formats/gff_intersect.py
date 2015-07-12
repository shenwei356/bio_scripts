#!/usr/bin/env python
# https://github.com/shenwei356/bio_scripts
# Author     : Wei Shen
# Contact    : shenwei356@gmail.com
# LastUpdate : 2015-06-26

from __future__ import print_function, division
import argparse
import os
import shutil
import sys
import gzip
from collections import defaultdict, Counter
from bx.intervals.intersection import Intersecter, Interval

parser = argparse.ArgumentParser(description="gff intersect",
                                 epilog="https://github.com/shenwei356/bio_scripts")

parser.add_argument('query', type=str, help='gff file b (query)')
parser.add_argument('subject', type=str, help='gff file a (subject)')
parser.add_argument('-e', '--embeded', action='store_true',
                    help='see what genes (subject) containing in specific regions (query)')
parser.add_argument('-s', '--split', action='store_true',
                    help='split results into multiple files')
parser.add_argument('-eu', '--extend-upstream', type=int, default=0,
                    help='extend N bases in the upstream [0]')
parser.add_argument('-ed', '--extend-downstream', type=int, default=0,
                    help='extend N bases in the downstream [0]')

args = parser.parse_args()

if args.extend_upstream and args.extend_upstream <= 0:
    sys.stderr.write('value of option --extend-upstream should be greater than 0\n')
    sys.exit(1)

if args.extend_downstream and args.extend_downstream <= 0:
    sys.stderr.write('value of option --extend-downstream should be greater than 0\n')
    sys.exit(1)

sys.stderr.write('building tree\n')
trees = dict()
with gzip.open(args.subject) if args.subject.endswith('.gz') else open(args.subject) as fh:
    genome = ''
    for line in fh:
        if line.isspace() or line[0] == '#':
            continue

        data = line.rstrip().split('\t')
        if len(data) != 9:
            sys.stderr.write('number of columns != 9: {}'.format(line))

        g, start, end, strand = data[0], int(data[3]), int(data[4]), data[6]
        if g != genome:
            genome = g
            trees[genome] = Intersecter()

        if strand == '+':
            start -= args.extend_upstream
            end += args.extend_downstream
        else:
            start -= args.extend_downstream
            end += args.extend_upstream

        if not args.embeded and strand == '-':  # complement strand
            start, end = -end, -start
        trees[genome].add_interval(Interval(start, end, value=data))

if args.split:
    outdir = '{}.intersect@{}'.format(os.path.normpath(os.path.basename(args.query)),
                                      os.path.normpath(os.path.basename(args.subject)))
    if os.path.exists(outdir):
        shutil.rmtree(outdir)
    os.makedirs(outdir)

sys.stderr.write('querying\n')
with gzip.open(args.query) if args.query.endswith('.gz') else open(args.query) as fh:
    for line in fh:
        if line.isspace() or line[0] == '#':
            continue
        data = line.rstrip().split('\t')
        if len(data) != 9:
            sys.stderr.write('number of columns != 9: {}'.format(line))

        genome, start, end, strand, product = data[0], int(data[3]), int(data[4]), data[6], data[8]

        if genome not in trees:
            continue

        overlaps = trees[genome].find(start, end)
        if len(overlaps) == 0:
            continue

        overlap_data, stats = list(), Counter()
        for x in overlaps:
            s, e = x.start, x.end
            if args.embeded:
                strand2 = '.'
            elif s > 0:
                strand2 = '+'
            else:  # complement strand
                s, e = -x.end, -x.start
                strand2 = '-'

            overlap, t = 0, ''
            if s <= start:
                if e >= end:
                    #   start ======== end
                    #     s ------------- e
                    if args.embeded:
                        continue
                    overlap = end - start + 1
                    t = 'cover'
                else:
                    #  start ======== end
                    #   s ------ e
                    if args.embeded:
                        continue
                    overlap = e - start + 1
                    t = 'overlap.upstream' if strand == '+' else 'overlap.downstream'
            else:
                if e >= end:
                    #   start ======== end
                    #           s ------ e
                    if args.embeded:
                        continue
                    overlap = end - s + 1
                    t = 'overlap.downstream' if strand == '+' else 'overlap.upstream'
                else:
                    #   start ======== end
                    #          s --- e
                    overlap = e - s + 1
                    t = 'embed'

            if args.embeded:
                frame = '.'
            elif strand == '+':
                frame = abs(s - start) % 3
            else:
                frame = abs(e - end) % 3

            stats[t] += 1
            if args.embeded:
                overlap_data.append(x.value)
            else:
                overlap_data.append([str(i) for i in
                                     [s, e, strand2, overlap, round(100 * overlap / (end - start + 1), 1), t, frame,
                                      x.value[-1]]])

        if args.split:
            fh_out = open(os.path.join(outdir, '{}_{}..{}..{} {}.gff'.format(genome, start, end, strand, product.replace('/','_'))),
                          'wt')
            fh_out.write('# {}'.format(line))
        else:
            fh_out = sys.stdout
            fh_out.write('>{}'.format(line))

        if args.embeded:
            sorted_overlap_data = sorted(overlap_data, key=lambda o: (o[0], o[1]))
        else:
            fh_out.write('# summary: {}\n'.format(stats))
            fh_out.write(
                '\t'.join(['start', 'end', 'strand', 'overlap', 'overlap%', 'type', 'frame', 'attribute']) + '\n')
            sorted_overlap_data = sorted(overlap_data, key=lambda o: (o[5], o[6], -float(o[4])))

        for overlap in sorted_overlap_data:
            fh_out.write('\t'.join(overlap) + '\n')
        fh_out.write('\n')

        if args.split:
            fh_out.close()
