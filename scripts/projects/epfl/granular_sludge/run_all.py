#!/usr/bin/env python2
# -*- coding: utf-8 -*-

"""
A script to run all the procedure on the granular sludge project.
"""

# Built-in modules #
import os
import shutil

# Internal modules #
import gefes

# First party modules #
from plumbing.processes import prll_map
from plumbing.timer     import Timer
from plumbing.autopaths import FilePath

# Third party modules #
from tqdm import tqdm

#################################### Load #####################################
# Load two projects #
proj1 = gefes.load("~/deploy/gefes/metadata/json/projects/epfl/granular_sludge_1/")
proj2 = gefes.load("~/deploy/gefes/metadata/json/projects/epfl/granular_sludge_2/")
projects = (proj1, proj2)
proj = proj1 + proj2

###############################################################################
print("# Get information for excel file #")
for s in proj: print s.short_name
for s in proj: print s.pair.fwd.count
for s in proj: print s.pair.rev.count
for s in proj: print s.pair.fwd.md5
for s in proj: print s.pair.rev.md5
for s in proj: print len(s.pair.fwd.first)

print("# Get special sample name #")
from collections import OrderedDict
d = open("/home/lucas/GEFES//raw/projects/epfl/granular_sludge/2/correspondances.csv").read()
d = OrderedDict((l.split()[0][:-13] + '.fastq.gz', l.split()[1]) for l in d.split('\n') if l)
for s in proj2: print s.short_name
for s in proj2: print d[s.pair.fwd.filename]