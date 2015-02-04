#!/usr/bin/env python3
# https://github.com/shenwei356/bio_scripts
# Author     : Wei Shen
# Contact    : shenwei356@gmail.com
# LastUpdate : 2015-02-04

import argparse
import csv
import logging
import sys
import re

# ===================================[ args ]=================================

parser = argparse.ArgumentParser(description="Join paired lines from two files into one file")

parser.add_argument("-v", "--verbose", help='Verbosely print information',
                    action="count", default=0)

parser.add_argument('infile1', type=argparse.FileType('r'),
                    help='Input file 1')
parser.add_argument('infile2', type=argparse.FileType('r'),
                    help='Input file 2')
parser.add_argument('outfile', nargs='*', type=argparse.FileType('w'),
                    default=sys.stdout, help='Output file')

parser.add_argument("-k", '--key', type=int, default=1,
                    help='Column number of key in csvfile')
parser.add_argument("-H", "--ignoretitle", help="Ignore title",
                    action="store_true")
parser.add_argument("-F", '--fs', type=str, default="\t",
                    help='Field separator [\\t]')
parser.add_argument("-Q", '--qc', type=str, default='"',
                    help='Quote char["]')

args = parser.parse_args()

# logging level
if args.verbose >= 2:
    logginglevel = logging.DEBUG
elif args.verbose == 1:
    logginglevel = logging.INFO
else:
    logginglevel = logging.WARN
logging.basicConfig(level=logginglevel,
                    format="[%(levelname)s] %(message)s")

logging.info("Column number of key in csvfile: {}".format(args.key))

# ===================================[ read csv ]=============================



def get_key_from_row(nrow, row):
    if nrow < args.key:
        logging.error(
            "-k ({}) is beyond number of column ({})".format(args.key, nrow))
        sys.exit(1)
    elif args.key < 1:
        args.key = 1
    key = row[args.key - 1].strip()
    return key


reader1 = csv.reader(iter(args.infile1.readline, ''), delimiter=args.fs, quotechar=args.qc)
reader2 = csv.reader(iter(args.infile2.readline, ''), delimiter=args.fs, quotechar=args.qc)

writer = csv.writer(args.outfile, delimiter=args.fs, quotechar=args.qc, quoting=csv.QUOTE_MINIMAL)

once = True
for row1, row2 in zip(reader1, reader2):
    if args.ignoretitle and once:  # Ignore title
        once = False
        continue

    nrow1, nrow2 = len(row1), len(row2)
    if nrow1 == 0 or nrow2 == 0:
        continue
    if nrow1 != nrow2:
        logging.error("unpaired column number: {} vs {}".format(nrow1, nrow2))
        sys.exit(1)

    key1, key2 = get_key_from_row(nrow1, row1), get_key_from_row(nrow2, row2)

    if key1 != key2:
        logging.error("keys do not match: {} vs {}".format(key1, key2))
        sys.exit(1)

    writer.writerow(row1)
    writer.writerow(row2)
