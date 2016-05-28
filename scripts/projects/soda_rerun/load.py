#!/usr/bin/env python2

"""
A script to load the main objects for the soda rerun project.
"""

# Built-in modules #
import os

# Internal modules #
import gefes

# First party modules #
from plumbing.processes import prll_map
from plumbing.timer     import Timer
from plumbing.cache     import LazyList

# Third party modules #
from tqdm import tqdm

# Constants #
user = os.environ.get('USER')

#################################### Load #####################################
# One project #
proj = gefes.projects['soda_rerun'].load()
samples = proj.samples
for s in samples: s.load()

# The bins #
bins = LazyList(lambda: proj.merged.results.binner.results.bins)
