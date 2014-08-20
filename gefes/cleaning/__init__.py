# Futures #
from __future__ import division

# Built-in modules #

# Internal modules #
from plumbing.autopaths import AutoPaths

# Third party modules #

###############################################################################
class Cleaner(object):
    """Takes care of cleaning the raw reads."""

    all_paths = """
    /fastqc/
    """

    def __repr__(self): return '<%s object of %s>' % (self.__class__.__name__, self.parent)
    def __len__(self): return len(self.pair)

    def __init__(self, pool):
        # Save parent #
        self.parent, self.pool = pool, pool
        # Auto paths #
        self.base_dir = self.parent.p.clean_dir
        self.p = AutoPaths(self.base_dir, self.all_paths)

    @property
    def ratio_discarded(self):
        return 1 - (len(self.pair) / len(self.pool.pair))