#!/usr/bin/env python2

"""
A script to run small snippets of code on the soda rerun project.
"""

# Modules #
import gefes
from plumbing.cache import LazyList

# One project #
proj = gefes.projects['soda_rerun'].load()
samples = proj.samples
for s in samples: s.load()

# The bins #
bins = LazyList(lambda: proj.merged.results.binner.results.bins)

###############################################################################
b = [b for b in bins if b.num == "1"][0]
b.report.generate()
print b.report.output_path