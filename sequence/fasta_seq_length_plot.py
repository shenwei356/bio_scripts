#!/usr/bin/env python
from __future__ import print_function

import sys
import os

from Bio import SeqIO

import seaborn as sns
import matplotlib as mpl
import matplotlib.pyplot as plt

usage = """
Usage: fasta_seq_length_plot.py fastafile [fastafile...]
"""

if len(sys.argv) <= 1:
    print(usage)
    sys.exit(0)

lengths = []

for file in sys.argv[1:]:
    if not os.path.exists(file):
        print("file not exists: %s" % file)
        sys.exit(0)

    with open(file + ".len", 'w') as fh:
        for seq in SeqIO.parse(file, "fasta"):
            l = len(seq)
            lengths.append(l)
            fh.write("%s\t%d\n" % (seq.id, l))

mpl.rc("figure", figsize=(8, 4))
sns.distplot(lengths)
plt.savefig(file + ".len.png")
