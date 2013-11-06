# Futures #
from __future__ import division

# Built-in modules #
from itertools import product

# Internal modules #
from gefes.common.autopaths import AutoPaths
from gefes.helper.clusterer import Clusterer
from gefes.common.cache import property_cached

# Third party modules #
import pandas

###############################################################################
class Binner(object):
    """Get a bunch of contigs, computes the coverage of the contigs,
    computes the tetranucleotide frequency, clusters them, and is
    responsible for delivering N bins, each with one or more contigs."""

    all_paths = """
    /mapping/
    """

    def __repr__(self): return '<%s object of %s>' % (self.__class__.__name__, self.parent)

    def __init__(self, aggregate):
        # Save parent #
        self.parent, self.aggregate = aggregate, aggregate
        # Auto paths #
        self.base_dir = self.parent.p.binning_dir
        self.p = AutoPaths(self.base_dir, self.all_paths)
        # Children #
        self.clusterer = Clusterer(self)
        # Output #
        self.bins = []

    @property_cached
    def frame(self):
        tetramers = ["".join(tetramer) for tetramer in product('ACGT', repeat=4)]
        columns = ['length'] + [s.id_name for s in self.aggregate] + ['freq_' +
        t for t in tetramers]

        frame = pandas.DataFrame([[c.length] +
            [s.mapper.coverage[c.name]["cov_mean"] for s in self.aggregate] +
            [c.tetra_nuc_freq.get(t, 0) for t in tetramers] for c
            in self.aggregate.assembly.contigs], columns=columns, index=[c.name
            for c in self.aggregate.assembly.contigs])

        return frame
