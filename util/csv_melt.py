#!/usr/bin/env python
# https://github.com/shenwei356/bio_scripts"
import argparse
import csv
import re
import sys
import pandas as pd

parser = argparse.ArgumentParser(
    description="Melt CSV file, you can append new column",
    epilog="https://github.com/shenwei356/bio_scripts")
parser.add_argument(
    'key',
    type=str,
    help=
    'Column name of key in csvfile. Multiple values shoud be separated by comma')
parser.add_argument('csvfile', type=str, help='CSV file with head row!')
parser.add_argument('--var_name',
                    type=str,
                    default='var_name',
                    help='name to use for the "variable" column')
parser.add_argument('--value_name',
                    type=str,
                    default='value_name',
                    help='name to use for the "value" column')
parser.add_argument('-a',
                    '--append',
                    type=str,
                    help='another column. format: column=value')
parser.add_argument('-o', '--outfile', type=str, help='output file')

parser.add_argument('--fs', type=str, default=",", help='field separator [,]')
parser.add_argument('--fs-out',
                    type=str,
                    help='field separator of ouput [same as --fs]')
parser.add_argument('--qc', type=str, default='"', help='Quote char["]')
parser.add_argument('-t',
                    action='store_true',
                    help='field separator is "\\t". Quote char is "\\t"')

args = parser.parse_args()

if args.t:
    args.fs, args.qc = '\t', '\t'
if not args.fs_out:
    args.fs_out = args.fs

pattern = '^([^=]+)=([^=]+)$'
if args.append:
    if not re.search(pattern, args.append):
        sys.stderr.write("bad format for option -a: {}".format(args.append))
        sys.exit(1)
    colname, colvalue = re.findall(pattern, args.append)[0]))

keys = list()
if ',' in args.key:
    keys = [k for k in args.key.split(',')]
else:
    keys = [args.key]

# ------------------------------------------------------------

df = pd.read_csv(args.csvfile,
                 sep=args.fs,
                 quotechar=args.qc)  # , index_col=keys)
df = pd.melt(df,
             id_vars=keys,
             var_name=args.var_name,
             value_name=args.value_name)
if args.append:
    df[colname] = pd.Series([colvalue] * len(df))

if args.outfile:
    df.to_csv(args.outfile, sep=args.fs, quotechar=args.qc, index=0)
else:
    df.to_csv(sys.stdout, sep=args.fs, quotechar=args.qc, index=0)
