# Futures #
from __future__ import division

# Built-in modules #

# Internal modules #
from gefes.common.autopaths import AutoPaths

# Third party modules #

###############################################################################
class Clusterer(object):
    """Recieves the matrix detailing all information for every contig,
    and is responsible for deciding which contigs go together."""

    all_paths = """
    /lorem.txt
    """

    def __repr__(self): return '<%s object of %s>' % (self.__class__.__name__, self.parent)

    def __init__(self, pool):
        # Save parent #
        self.parent, self.pool = pool, pool
        # Auto paths #
        self.base_dir = self.parent.p.clustering
        self.p = AutoPaths(self.base_dir, self.all_paths)

    def run(self):
        pass