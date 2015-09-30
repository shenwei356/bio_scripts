#!/usr/bin/env python
# -*- coding: utf-8 -*-
# https://github.com/shenwei356/bio_scripts
from __future__ import print_function
import argparse
import logging
import os
import sys
import gzip
import re
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def parse_args():
    parser = argparse.ArgumentParser(description="Translate DNA to peptide")

    parser.add_argument("-v", "--verbose", help='verbosely print information',
                        action="count", default=0)

    group = parser.add_mutually_exclusive_group()
    group.add_argument("--stdin", action="store_true",
                       help='read from stdin, one sequence per line')
    group.add_argument('-i', '--infile', type=str,
                       help='file name should like this: infile.[fasta|fa|fastq|fq][.gz]')
    parser.add_argument('-f', '--format', type=str, # default='fasta',
                        help='seqence format: fasta |fastq  [fasta]')
    parser.add_argument('-t', '--table', type=int, default=1,
                        help='genetic code table (detail: http://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi ) [1]')

    args = parser.parse_args()
    if not (args.stdin or args.infile):
        sys.stderr.write("option --stdin or -i should be given\n")
        sys.exit(1)
    if args.format and not args.format in ['fasta', 'fastq']:
        sys.stderr.write("option -f | --format should be 'fasta' or 'fastq'\n")
        sys.exit(1)
    if args.stdin and not args.format:
        sys.stderr.write("option -f | --format should be given when --stdin is set.\n")
        sys.exit(1)

    return args


if __name__ == '__main__':
    args = parse_args()

    file, seq_format, fh = args.infile, args.format,  None,
    if file:
        if not seq_format:
            found = re.search(r'(?i)(fasta|fa|fastq|fq)(.gz)?$', file)
            if not found:
                print("invalid file name suffix.\nfile name should like this: infile.[fasfa|fa|fastq|fq][.gz]",
                      file=sys.stderr)
                sys.exit(1)
            seq_format, is_gz = found.groups()
            if seq_format == 'fa':
                seq_format = 'fasta'
            if seq_format == 'fq':
                seq_format = 'fastq'

        fh = gzip.open(file, 'rt') if file.endswith('.gz') else open(file, 'r')
    else:
        fh = sys.stdin
        seq_format = args.format

    for seq in SeqIO.parse(fh, seq_format):
        SeqIO.write([SeqRecord(seq.seq.translate(table=args.table), id=seq.id, description=seq.description)], sys.stdout, 'fasta')

    fh.close()
