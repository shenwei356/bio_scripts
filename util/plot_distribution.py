#!/usr/bin/env python
from __future__ import print_function

import sys
import os


import seaborn as sns
import matplotlib as mpl
import matplotlib.pyplot as plt


mpl.rc("figure", figsize=(8, 4))
sns.distplot(lengths)
plt.savefig(file + ".len.png")
