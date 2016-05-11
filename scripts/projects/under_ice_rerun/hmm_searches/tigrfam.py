#!/usr/bin/env python2

"""
A script to generate all the reports of the under_ice_rerun project.
"""

# Get the deploy version #
import os, sys
home = os.environ.get('HOME', '~') + '/'
sys.path.insert(0, home + 'deploy/gefes')

# Modules #
import gefes
from tqdm import tqdm

# Three projects #
bt = gefes.projects['under_ice_rerun_bt'].load()
lb = gefes.projects['under_ice_rerun_lb'].load()
kt = gefes.projects['under_ice_rerun_kt'].load()
for s in bt.samples: s.load()
for s in lb.samples: s.load()
for s in kt.samples: s.load()
samples = tuple(bt.samples + lb.samples + kt.samples)
projects = (bt, lb, kt)

###############################################################################
# Thread pool #
from multiprocessing import Pool
def hmmsearch(bin): return bin.tigrfams.run(cpus=1)
bins = kt.merged.results.binner.results.bins + \
       lb.merged.results.binner.results.bins + \
       bt.merged.results.binner.results.bins
pool = Pool(processes=64)
iterator = pool.imap(hmmsearch, bins, chunksize=1)
i = 0
for result in tqdm(iterator, total=len(bins)):
    print "Done bin %s." % bins[i]
    i += 1