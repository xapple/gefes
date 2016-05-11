#!/usr/bin/env python2

"""
A script to generate all the reports of the under_ice_rerun project.
"""

# Modules #
import gefes, shutil
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
# Merged #
for p in projects: print p.merged.report.generate()

# Samples #
print "Generating reports for %i samples" % len(samples)
for s in tqdm(samples):
    s.report.generate()

# Basic assemblies #
assemblies = [(a,p) for p in projects for a in p.assemblies.values()]
for a,p in tqdm(assemblies):
    print "Report for project '%s', assembly '%s'" % (p.name, a)
    a.report.generate()

# Copy all #
assemblies = [(a,p) for p in projects for a in p.assemblies.values()]
for a,p in assemblies: shutil.copy(a.report.output_path, a.report.copy_base)
for s in samples:      shutil.copy(s.report.output_path, s.report.copy_base)
