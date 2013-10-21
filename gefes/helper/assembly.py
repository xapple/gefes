# Futures #
from __future__ import division

# Built-in modules #

# Internal modules #
from gefes.common.autopaths import AutoPaths
from gefes.graphs import assembly_plots

# Third party modules #

###############################################################################
class Assembly(object):
    """An assembly analysis."""

    all_paths = """
    /graphs/
    /velvet/contigs.fasta
    /amos/contigs.afg
    """

    def __repr__(self): return '<%s object of %s>' % (self.__class__.__name__, self.parent)

    def __init__(self, pool):
        # Save parent #
        self.parent, self.pool = pool, pool
        # Auto paths #
        self.base_dir = self.parent.base_dir
        self.p = AutoPaths(self.base_dir, self.all_paths)

    def assemble(self):
        pass

    def make_plots(self):
        for cls_name in assembly_plots.__all__:
            cls = getattr(assembly_plots, cls_name)
            cls(self).plot()