#!/usr/bin/python
# -*- coding: utf-8 -*-

import sys
import os

if len(sys.argv) != 2:
    print "\nUsage: %s  <alignment file>\n" % os.path.basename(sys.argv[0])
    sys.exit(1)

blast          = sys.argv[1]

with open(blast, 'r') as fp:
    init = ""
    for line in fp:
        if not line.startswith("#"):
            item         = line.strip().split("\t")
            if init != item[0]:
                print line.strip()
                init = item[0]
