#!/usr/bin/env python2

"""
A script to generate all the reports of the under_ice_rerun project.
"""

# Modules #
import gefes, tqdm

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
print "Generating report on %i samples" % len(samples)
for s in tqdm(samples):
    s.report.generate()

for p in projects:
    for a in p.assemblies.values():
        print "Report for project '%s', assembly '%s'" % (p.name, a)
        a.report.generate()