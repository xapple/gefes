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
# One contig #
#for c in tqdm(bt.merged.results.contigs): break
#c = bt.merged.results.contigs[0]
#print c.proteins.results.faa.count
#c.pfams
#from plumbing.timer import Timer
#with Timer(): c.pfams.run(cpus=4)
##Start at: 2016-05-11 17:35:28.055135
##End at: 2016-05-11 17:37:56.227389
##Total elapsed time: 0:02:28.172431
#
## One bin #
#for b in tqdm(bt.merged.results.binner.results.bins): break
#b = bt.merged.results.binner.results.bins[0]
#print b.faa.count
#b.pfams
#from plumbing.timer import Timer
#with Timer(): b.pfams.run(cpus=1)
##Start at: 2016-05-11 17:49:48.018523
##End at: 2016-05-11 17:54:08.518643
##Total elapsed time: 0:04:20.500248

# Thread pool #
from multiprocessing import Pool
def hmmsearch(bin): return bin.pfams.run(cpus=1)
bins = kt.merged.results.binner.results.bins + \
       lb.merged.results.binner.results.bins + \
       bt.merged.results.binner.results.bins
pool = Pool(processes=50)
iterator = pool.imap(hmmsearch, bins, chunksize=1)
i = 0
for result in tqdm(iterator, total=len(bins)):
    print "Done bin %s." % bins[i]
    i += 1