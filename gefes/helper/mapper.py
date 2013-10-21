# Futures #
from __future__ import division

# Built-in modules #

# Internal modules #
from gefes.common import AutoPaths

# Third party modules #

###############################################################################
class Mapper(object):
    """Lorem ipsum."""

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
