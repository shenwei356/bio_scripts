#!/usr/bin/env python
from __future__ import print_function

import argparse
import sys
import os

import seaborn as sns
import matplotlib as mpl
import matplotlib.pyplot as plt


parser = argparse.ArgumentParser(description='Plot distribution',
                                 epilog="https://github.com/shenwei356/bio_scripts")

parser.add_argument('-i', '--infile', nargs='?', type=argparse.FileType('r'),
                    default=sys.stdin,  help='Input file')

parser.add_argument('-o', '--outfile', nargs='?', type=str,
                    default='dist.png', help='Output file')

parser.add_argument('--width', type=int, default=8, help='Figure width')
parser.add_argument('--height', type=int, default=6, help='Figure heigth')

parser.add_argument(
    '-t', '--title', type=str, default='Distribution Plot', help='Figure Title')
parser.add_argument(
    '-x', '--xlabel', type=str, default='Value', help='Figure X label')
parser.add_argument(
    '-y', '--ylabel', type=str, default='Frequency', help='Figure Y label')

args = parser.parse_args()

data = []
for line in args.infile:
    data.append(float(line.strip()))

mpl.rc("figure", figsize=(args.width, args.height))

figure = sns.distplot(data)

figure.set_title(args.title)
figure.set_xlabel(args.xlabel)
figure.set_ylabel(args.ylabel)

plt.savefig(args.outfile)
