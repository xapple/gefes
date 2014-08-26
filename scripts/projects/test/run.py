#!/usr/bin/env python2

"""
A script to contain the procedure for running the test sample.
"""

# Built-in modules #

# Internal modules #
import gefes
from gefes.groups.aggregate import Aggregate

# Third party modules #

# Constants #
proj = gefes.projects['test']
samples = proj.samples

# Global settings #
proj.kmer_size = 81

################################ Preprocessing ################################
for s in samples:
    s.pair.fwd.fastqc.run()
    s.pair.rev.fastqc.run()
    s.quality_checker.run()
    s.clean.fwd.fastqc.run()
    s.clean.rev.fastqc.run()
    s.clean.fwd.graphs['LengthDist'].plot()
    s.clean.rev.graphs['LengthDist'].plot()
    s.pair.fwd.avg_quality
    s.pair.rev.avg_quality
    s.report.generate()

################################### Aggregate ##################################
hypolimnion = Aggretate()
metalimnion =
epilimnion =


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