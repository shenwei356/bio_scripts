#!/usr/bin/env python3
# https://github.com/shenwei356/bio_scripts

import argparse
import sys
from collections import defaultdict, Counter
import pysam

parser = argparse.ArgumentParser(description="bam2gff",
                                 epilog="https://github.com/shenwei356/bio_scripts")

parser.add_argument('bamfile', type=str, help='bam file')
parser.add_argument('-c', '--cache-size', type=int, default=1000, help='cache size [1000]')
parser.add_argument('-m', '--match-proportion', type=float, default=0.75,
                    help='minimum match proportion to define properly paired ends [0.75]')

parser.add_argument("-v", "--verbose", help='verbosely print information',
                    action="count", default=0)

args = parser.parse_args()

pairs = defaultdict(lambda: defaultdict(dict))
stats = Counter()
samfile = pysam.AlignmentFile(args.bamfile, "rb")
for read in samfile.fetch():
    if read.is_proper_pair and not read.is_secondary:
        if read.reference_length < read.query_length * 0.75:  # full match
            stats['bad match'] += 1
            continue
        key = '_'.join([str(x) for x in sorted([read.reference_start, read.next_reference_start])])
        pairs[read.query_name][key]['read1' if read.is_read1 else 'read2'] = {'start': read.reference_start,
                                                                              'end': read.reference_end,
                                                                              'ref': samfile.getrname(
                                                                                  read.reference_id),
                                                                              'reverse': read.is_reverse}

        if 'read1' in pairs[read.query_name][key] and 'read2' in pairs[read.query_name][key]:
            read1, read2 = pairs[read.query_name][key]['read1'], pairs[read.query_name][key]['read2']

            if not read1['reverse']:
                strand, start, end = '+', read1['start'], read2['end']
            else:
                strand, start, end = '-', read2['start'], read1['end']

            sys.stdout.write('\t'.join(
                [read1['ref'], 'bam2gff.py', 'paired_ends', str(start + 1), str(end), '.', strand, '.',
                 read.query_name]) + "\n")

            stats['paired'] += 1

            del pairs[read.query_name][key]

samfile.close()

for query, sites in pairs.items():
    if len(sites) == 0:
        continue
    stats['unpaired'] += 1

sys.stderr.write('{} summary: {}\n'.format(args.bamfile, stats))
