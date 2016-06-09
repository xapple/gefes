#!/usr/bin/env python2

"""
A script to run small snippets of code on the under ice rerun project.
"""

import names

# Modules #
import gefes, shutil

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
for proj in projects: proj.merged.results.trait_annotations.run()
