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
from Bio import motifs
from Bio import SeqIO
from Bio.Seq import Seq


def parse_args():
    parser = argparse.ArgumentParser(description="Position Weight Matrices of sequence")

    parser.add_argument("-v", "--verbose", help='verbosely print information',
                        action="count", default=0)

    group = parser.add_mutually_exclusive_group()
    group.add_argument("--stdin", action="store_true",
                       help='read from stdin, one sequence per line')
    group.add_argument('-i', '--infile', type=str,
                       help='file name should like this: infile.[fasta|fa|fastq|fq][.gz]')

    args = parser.parse_args()
    if not ( args.stdin or args.infile ):
        sys.stderr.write("option --stdin or -i should be given\n")
        sys.exit(1)

    return args


def seq_iter(file):
    if file:
        found = re.search(r'(?i)(fasta|fa|fastq|fq)(.gz)?$', file)
        if not found:
            sys.stderr.write("invalid file name suffix.\nfile name should like this: infile.[fasfa|fa|fastq|fq][.gz]\n")
            sys.exit(1)
        seq_format, is_gz = found.groups()
        if seq_format == 'fa':
            seq_format = 'fasta'
        if seq_format == 'fq':
            seq_format = 'fastq'

        fh = gzip.open(file, 'rt') if is_gz else open(file, 'r')
        for record in SeqIO.parse(fh, seq_format):
            yield record.seq
        fh.close()
    else:
        for line in sys.stdin:
            yield Seq(line.strip())


if __name__ == '__main__':
    args = parse_args()
    seqs = seq_iter(args.infile)
    seqs2 = [seq for seq in seqs if not 'N' in seq]
    m = motifs.create(seqs2)
    print(m.pwm)
    # print(m.pssm)
    # m.weblogo("motif.png")
