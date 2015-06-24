#!/usr/bin/env python3
# https://github.com/shenwei356/bio_scripts

import argparse
import sys
from collections import defaultdict
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

pairs = defaultdict(dict)


def pairs2gff(pairs):
    todelete = list()
    for query, pair in pairs.items():
        if len(pair.keys()) < 2:
            continue
        read1, read2 = pair['read1'], pair['read2']

        if not read1['reverse']:
            strand, start, end = '+', read1['start'], read2['end']
        else:
            strand, start, end = '-', read2['start'], read1['end']

        sys.stdout.write('\t'.join(
            [query, 'bam2gff.py', 'paired_ends', str(start + 1), str(end), '.', strand, '.',
             read1['ref']]) + "\n")
        todelete.append(query)
    for query in todelete:
        del pairs[query]


samfile = pysam.AlignmentFile(args.bamfile, "rb")
for read in samfile.fetch():
    if read.is_proper_pair and not read.is_secondary:
        if read.reference_length < read.query_length * 0.75:  # full match
            continue

        pairs[read.query_name]['read1' if read.is_read1 else 'read2'] = {'start': read.reference_start,
                                                                         'end': read.reference_end,
                                                                         'ref': samfile.getrname(read.reference_id),
                                                                         'reverse': read.is_reverse}
        if len(pairs) > args.cache_size:
            pairs2gff(pairs)

samfile.close()

pairs2gff(pairs)
