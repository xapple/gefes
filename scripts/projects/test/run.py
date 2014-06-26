#!/usr/bin/env python2

"""
A script to contain examples commands for running the pipeline.
"""

# Don't run it #
import sys
sys.exit("Copy paste the commands you want in ipython, don't run this script.")

# Modules #
import gefes

# Constants #
proj = gefes.projects['test']
pools = proj.pools

################################### Cleaning ##################################
self.joiner = Joiner(self)
# Check quality #
self.quality_checker = QualityChecker(self)
# Final files #
self.fwd = self.quality_checker.fwd
self.rev = self.quality_checker.rev
self.pair = self.quality_checker.paired

def clean():
    # Fastqc on the pools #
    for p in pools: p.pair.fastqc(directory=p.p.fastqc_dir)
    # Clean the pools #
    for p in pools: p.clean_reads()
    # The clean graphs #
    for p in pools: p.make_plots()
    # Fastqc on the result #
    for p in pools: p.cleaner.pair.fastqc(directory=p.cleaner.p.fastqc_dir)

################################### Cleaning ##################################
# Assemble locally #
proj.assemble()
# The assembly graphs #
proj.assembly.graphs[0].plot()
# Index the result #
proj.assembly.index()
# Map the reads #
for p in pools: p.mapper.map()

# Binning frame #
proj.binner.export_frame()
# Clustering #
proj.binner.clusterer.run()