#!/usr/bin/python
# from https://github.com/jameslz/benchmark
import sys

if len(sys.argv) != 4:
    print "\nUsage: %s  <alignment> <evalue> <bit_score>\n" % sys.argv[0]
    sys.exit(1)

blast          = sys.argv[1]
max_evalue     = float(sys.argv[2])
min_bit_score  = float(sys.argv[3])

with open(blast, 'r') as fp:
    init = ""
    for line in fp:
        if not line.startswith("#"):
            item         = line.strip().split("\t")
            evalue       = float(item[10])
            bit_score    = float(item[11])
            if init != item[0]:
                if evalue <= max_evalue and bit_score >= min_bit_score:
                    print line.strip()
                init = item[0] 
