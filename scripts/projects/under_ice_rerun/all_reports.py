#!/usr/bin/env python2

"""
A script to generate all the reports of the under_ice_rerun project.
"""

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
# Merged #
for p in projects:
    print "Generating report for merged assembly %s" % p
    print p.merged.report.generate()

# Samples #
print "Generating reports for %i samples" % len(samples)
for s in tqdm(samples):
    s.report.generate()

# Basic assemblies #
assemblies = [(a,p) for p in projects for a in p.assemblies.values()]
print "Generating reports for %i assemblies" % len(assemblies)
for a,p in tqdm(assemblies):
    #print "Report for project '%s', assembly '%s'" % (p.name, a)
    a.report.generate()

# Samples clean cache #
if False:
    for s in samples: s.report.cache_dir.remove()
    for s in samples: s.report.cache_dir.create()

# Sample clean graphs #
for s in samples:
    s.mono_mapper.results.graphs.mean_coverage.path.remove()
    s.mono_mapper.results.graphs.percent_covered.path.remove()
    s.mapper_merged.results.graphs.mean_coverage.path.remove()
    s.mapper_merged.results.graphs.percent_covered.path.remove()