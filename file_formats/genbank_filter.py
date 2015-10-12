#!/usr/bin/env python
from __future__ import print_function
import re
import sys
import argparse

parser = argparse.ArgumentParser(
    description='filter gene records by regular expression from genbank file',
    epilog="https://github.com/shenwei356/bio_scripts")
parser.add_argument('gbfile', type=str, help='genbank file')
parser.add_argument('pattern',
                    type=str,
                    help='pattern (regular expression) [.]')
args = parser.parse_args()

with open(args.gbfile) as fh:
    tmp = ''
    for line in fh:
        if line.startswith('     gene '):
            if tmp == '':
                tmp = line
            else:
                if re.search(args.pattern, tmp):
                    sys.stdout.write(tmp)
                tmp = line
        elif line != '':
            tmp += line
if re.search(args.pattern, tmp):
    sys.stdout.write(tmp)
