# Futures #
from __future__ import division

# Built-in modules #

# Internal modules #
from gefes.common import AutoPaths
from gefes.helper.mapper import Mapper
from gefes.helper.linkage import Linkage
from gefes.helper.clusterer import Clusterer

# Third party modules #

###############################################################################
class Binner(object):
    """Get a bunch of contigs, computes the coverage of the contigs,
    computes the tetranucleotide frequency, clusters them, and is
    responsible for delivering N bins, each with one or more contigs."""

    all_paths = """
    /lorem.txt
    """

    def __repr__(self): return '<%s object of %s>' % (self.__class__.__name__, self.parent)

    def __init__(self, pool):
        # Save parent #
        self.parent, self.pool = pool, pool
        # Auto paths #
        self.base_dir = self.parent.base_dir
        self.p = AutoPaths(self.base_dir, self.all_paths)
        # Children #
        self.mapper = Mapper(self)
        self.linkage = Linkage(self)
        self.clusterer = Clusterer(self)
        # Output #
        self.bins = []

    def bam_to_coverage(self):
        pass

    def add_nuc_freq_info(self):
        pass