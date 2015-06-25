#!/usr/bin/env python
# https://github.com/shenwei356/bio_scripts

from __future__ import print_function
import argparse
import sys
import gzip
from collections import defaultdict
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

parser = argparse.ArgumentParser(description="extract_cds_by_gff")
parser.add_argument('-t', '--type', type=str,
                    default='CDS', help='gene type. "." for any types. [CDS]')
parser.add_argument('-us', '--up-stream', type=int,
                    default=0, help='up stream length [0]')
parser.add_argument('-ds', '--down-stream', type=int,
                    default=0, help='down stream length [0]')
parser.add_argument('gff_file', type=str, help='gff file')
parser.add_argument('fasta_file', type=str, help='fasta file')
args = parser.parse_args()
if not (args.up_stream >= 0 and args.down_stream >= 0):
    print('value of --up-stream and --down-stream should be >= 0', file=sys.stderr)
    sys.exit(1)


def read_gff_file(file):
    genes = defaultdict(list)
    with open(file, 'rt') as fh:
        for row in fh:
            data = row.strip().split('\t')
            if len(data) < 9:
                continue
            name = data[0]
            gene = dict()
            gene['type'], gene['start'], gene['end'], gene['strand'] = data[2], int(data[3]), int(data[4]), data[6]
            genes[name].append(gene)

    return genes


genes = read_gff_file(args.gff_file)

fh = gzip.open(args.fasta_file, 'rt') if args.fasta_file.endswith('.gz') else open(args.fasta_file, 'r')
for record in SeqIO.parse(fh, 'fasta'):
    name, genome = record.id, record.seq

    if name not in genes:
        continue

    for gene in genes[name]:
        if args.type != '.' and gene['type'].lower() != args.type.lower():
            continue
        seq = ''
        if gene['strand'] == '+':
            s = gene['start'] - args.up_stream - 1
            s = 0 if s < 0 else s
            seq = genome[s:gene['end'] + args.down_stream]
        else:
            s = gene['start'] - args.down_stream - 1
            s = 0 if s < 0 else s
            seq = genome[s:gene['end'] + args.up_stream].reverse_complement()
        SeqIO.write(
            SeqRecord(seq, id='{}_{}..{}..{}'.format(name, gene['start'], gene['end'], gene['strand']), description=''),
            sys.stdout, 'fasta')
fh.close()
