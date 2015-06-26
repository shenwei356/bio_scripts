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

parser = argparse.ArgumentParser(description="Split paired lines into two files")

parser.add_argument('csvfile', nargs='*', type=argparse.FileType('r'),
                    default=sys.stdin, help='Input file(s)')
parser.add_argument("-v", "--verbose", help='Verbosely print information',
                    action="count", default=0)

parser.add_argument('outfile1', type=argparse.FileType('w'),
                    default="out_1.tab", help='Output file 1')
parser.add_argument('outfile2', type=argparse.FileType('w'),
                    default="out_2.tab", help='Output file 2')

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

writer1 = csv.writer(args.outfile1, delimiter=args.fs, quotechar=args.qc, quoting=csv.QUOTE_MINIMAL)
writer2 = csv.writer(args.outfile2, delimiter=args.fs, quotechar=args.qc, quoting=csv.QUOTE_MINIMAL)

cnt, sum = 0, 0
stdinflag = False

# If "iter(sys.stdin.readline, '')" in the flowing for-loop, first line
# of stdin will be missing
if args.csvfile is sys.stdin:
    logging.info("read data from STDIN")
    stdinflag = True
    args.csvfile = [iter(sys.stdin.readline, '')]


def get_key_from_row(nrow, row):
    if nrow < args.key:
        logging.error(
            "-k ({}) is beyond number of column ({})".format(args.key, nrow))
        sys.exit(1)
    elif args.key < 1:
        args.key = 1
    key = row[args.key - 1].strip()
    return key


key0, row0, flag = '', '', True

for f in args.csvfile:
    if not stdinflag:
        logging.info("read data from file")
        f = iter(f.readline, '')
    reader = csv.reader(f, delimiter=args.fs, quotechar=args.qc)

    once = True
    for row in reader:
        if args.ignoretitle and once:  # Ignore title
            once = False
            continue

        nrow = len(row)
        if nrow == 0:
            continue

        sum += 1
        key = get_key_from_row(nrow, row)

        if key0 == '':
            key0, row0 = key, row
            continue

        if flag:
            if key0 != key:
                logging.error("unpaired key: line {} {} vs line {} {} ".format(sum - 1, row0, sum, row))
                sys.exit(1)
            else:
                writer1.writerow(row0)
                writer2.writerow(row)

        flag = not flag
        key0, row0 = key, row

if flag:
    logging.error("unpaired record remain: {}".format(row0))
    sys.exit(1)