#!/usr/bin/env python
from __future__ import print_function

import sys
import os

from Bio import SeqIO
from Bio.SeqUtils import GC

import seaborn as sns
import matplotlib as mpl
import matplotlib.pyplot as plt

usage = """
Usage: fasta_seq_gc_content_plot.py fastafile [fastafile...]
"""

if len(sys.argv) <= 1:
    print(usage)
    sys.exit(0)

gc = []

for file in sys.argv[1:]:
    if not os.path.exists(file):
        print("file not exists: %s" % file)
        sys.exit(0)

    with open(file + ".gc", 'w') as fh:
        for seq in SeqIO.parse(file, "fasta"):
            gccontent = GC(seq.seq)
            gc.append(gccontent)
            fh.write("%s\t%d\n" % (seq.id, gccontent))

mpl.rc("figure", figsize=(8, 4))
sns.distplot(gc)
plt.savefig(file + ".gc.png")
