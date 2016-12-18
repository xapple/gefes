#!/usr/bin/env python2

"""
A script to load the main objects for the under_ice rerun project.

ipython -i ~/repos/gefes/scripts/projects/under_ice_rerun/load.py
"""

# Built-in modules #
import os, shutil

# Internal modules #
import gefes
from gefes.groups.lump import Lump

# First party modules #
from plumbing.processes import prll_map
from plumbing.timer     import Timer
from plumbing.cache     import LazyList

# Third party modules #
from tqdm import tqdm

# Constants #
user = os.environ.get('USER')

#################################### Load #####################################
# Three projects #
bt = gefes.projects['under_ice_rerun_bt'].load()
lb = gefes.projects['under_ice_rerun_lb'].load()
kt = gefes.projects['under_ice_rerun_kt'].load()
for s in bt.samples: s.load()
for s in lb.samples: s.load()
for s in kt.samples: s.load()
samples = tuple(bt.samples + lb.samples + kt.samples)
projects = (bt, lb, kt)

# Special lump #
lump = Lump("ice_lump", projects)

# The bins #
bins = LazyList(lambda: kt.merged.results.binner.results.bins + \
                        lb.merged.results.binner.results.bins + \
                        bt.merged.results.binner.results.bins)
