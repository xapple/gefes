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
    pool.pair.fastqc(out_dir=pool.p.fastqc_dir)
    pool.quality_checker = QualityChecker(pool.pair, pool.clean)
    pool.quality_checker.run()
    pool.clean.fastqc(out_dir=pool.p.fastqc_dir)

################################### Assembly ##################################
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