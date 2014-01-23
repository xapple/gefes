# Futures #
from __future__ import division

# Built-in modules #

# Internal modules #
from gefes.common.autopaths import AutoPaths
from gefes.helper.clusterer import Clusterer
from gefes.common.cache import property_cached
from gefes.fasta.single import FASTA
from gefes.helper.contig import Contig


# Third party modules #
###############################################################################
class Binner(object):
    """Makes the matrix detailing all information for every contig,
    The output is N bins, each with one or more contigs inside."""

    all_paths = """
    /clustering/
    /bins/
    /frame.csv
    """

    def __repr__(self): return '<%s object of %s>' % (self.__class__.__name__, self.parent)

    def __init__(self, aggregate):
        # Save parent #
        self.parent, self.aggregate = aggregate, aggregate
        # Auto paths #
        self.base_dir = self.parent.p.binning_dir
        self.p = AutoPaths(self.base_dir, self.all_paths)
        # Children #
        self.clusterer=Clusterer(self)
        # Output #
        self.bins = []


    def run():
        if self.clusterer.clusters is None: self.clusterer.run()
