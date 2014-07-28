#!/usr/bin/env python2

"""
A script to contain the procedure for running the test sample.
"""

# Built-in modules #

# Internal modules #
import gefes
from gefes.cleaning.quality import QualityChecker

# Third party modules #

# Constants #
proj = gefes.projects['test']
pools = proj.pools

# Global settings #
proj.kmer_size = 81

################################### Cleaning ##################################
for pool in pools:
    # Check quality #
    pool.quality_checker = QualityChecker(pool)
    # Final files #
    pool.fwd  = pool.quality_checker.fwd
    pool.rev  = pool.quality_checker.rev
    pool.pair = pool.quality_checker.paired

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