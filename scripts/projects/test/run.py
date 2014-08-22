#!/usr/bin/env python2

"""
A script to contain the procedure for running the test sample.
"""

# Built-in modules #

# Internal modules #
import gefes
from gefes.cleaning.quality import QualityChecker
from gefes.report.sample import SampleReport

# Third party modules #

# Constants #
proj = gefes.projects['test']
samples = proj.samples

# Global settings #
proj.kmer_size = 81

################################ Preprocessing ################################
for s in samples:
    s.pair.fastqc(out_dir=s.p.fastqc_dir)
    s.quality_checker = QualityChecker(s.pair, s.clean)
    s.quality_checker.run()
    s.clean.fastqc(out_dir=s.p.fastqc_dir)

#################################### Report ###################################
for s in samples:
    s.report = SampleReport(s)
    s.report.generate()

################################### Assembly ##################################
# Assemble locally #
proj.assemble()
# The assembly graphs #
proj.assembly.graphs[0].plot()
# Index the result #
proj.assembly.index()
# Map the reads #
for s in samples: s.mapper.map()

# Binning frame #
proj.binner.export_frame()
# Clustering #
proj.binner.clusterer.run()