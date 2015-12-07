#!/usr/bin/env python
# -*- coding: utf-8 -*-
# https://github.com/shenwei356/bio_scripts

from __future__ import print_function

import argparse
import gzip
import sys
from collections import defaultdict

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

parser = argparse.ArgumentParser(description="extract_cds_by_gff")
parser.add_argument('-t',
                    '--type',
                    type=str,
                    default='CDS',
                    help='gene type. "." for any types. [CDS]')
parser.add_argument('-us',
                    '--up-stream',
                    type=int,
                    default=0,
                    help='up stream length [0]')
parser.add_argument('-ds',
                    '--down-stream',
                    type=int,
                    default=0,
                    help='down stream length [0]')
parser.add_argument('-j',
                    '--just',
                    action="store_true",
                    help='only output up and down stream')
parser.add_argument('gff_file', type=str, help='gff file')
parser.add_argument('fasta_file', type=str, help='fasta file')
args = parser.parse_args()
if not (args.up_stream >= 0 and args.down_stream >= 0):
    print('value of --up-stream and --down-stream should be >= 0',
          file=sys.stderr)
    sys.exit(1)
if args.just:
    if args.up_stream and args.down_stream or not (args.up_stream or
                                                       args.down_stream):
        print(
            'when using option --just, ONE of --up-stream and --down-stream should given',
            file=sys.stderr)
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
            gene['type'], gene['start'], gene['end'], gene['strand'], gene[
                'product'
            ] = data[2], int(data[3]), int(
                data[4]), data[6], data[8]
            genes[name].append(gene)

    return genes


genes = read_gff_file(args.gff_file)

fh = gzip.open(args.fasta_file,
               'rt') if args.fasta_file.endswith('.gz') else open(
                   args.fasta_file, 'r')
for record in SeqIO.parse(fh, 'fasta'):
    name, genome = record.id, record.seq
    genomesize = len(genome)
    if name not in genes:
        continue

    for gene in genes[name]:
        if args.type != '.' and gene['type'].lower() != args.type.lower():
            continue
        seq = ''
        flag = ''
        if gene['strand'] == '+':
            if args.just:
                if args.up_stream:
                    s = gene['start'] - args.up_stream - 1
                    e = gene['start'] - 1
                    flag = 'jus..{}'.format(args.up_stream)
                else:
                    s = gene['end']
                    e = gene['end'] + args.down_stream
                    flag = 'jds..{}'.format(args.down_stream)
            else:
                s = gene['start'] - args.up_stream - 1
                s = 0 if s < 0 else s
                e = gene['end'] + args.down_stream
                if args.up_stream:
                    flag = 'us..{}'.format(args.up_stream)
                else:
                    flag = 'ds..{}'.format(args.down_stream)

            s = 0 if s < 0 else s
            end = genomesize - 1 if e > genomesize - 1 else e
            seq = genome[s:e]
        else:
            if args.just:
                if args.up_stream:
                    s = gene['end']
                    e = gene['end'] + args.up_stream
                    flag = 'jus..{}'.format(args.up_stream)
                else:
                    s = gene['start'] - args.down_stream - 1
                    e = gene['start'] - 1
                    flag = 'jds..{}'.format(args.down_stream)
            else:
                s = gene['start'] - args.down_stream - 1
                s = 0 if s < 0 else s
                e = gene['end'] + args.up_stream
                if args.up_stream:
                    flag = 'us..{}'.format(args.up_stream)
                else:
                    flag = 'ds..{}'.format(args.down_stream)

            s = 0 if s < 0 else s
            end = genomesize - 1 if e > genomesize - 1 else e
            seq = genome[s:e].reverse_complement()

        if args.up_stream or args.down_stream:
            id = '{}_{}..{}..{}_{}'.format(name, gene['start'], gene['end'],
                                           gene['strand'], flag)
        else:
            id = '{}_{}..{}..{}'.format(name, gene['start'], gene['end'],
                                        gene['strand'])
        SeqIO.write(
            SeqRecord(seq,
                      id=id,
                      description=gene['product']),
            sys.stdout,
            'fasta')
fh.close()
