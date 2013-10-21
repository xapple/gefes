# Futures #
from __future__ import division

# Built-in modules #

# Internal modules #
from gefes.common import flatten
from gefes.common.autopaths import AutoPaths
from gefes.graphs import assembly_plots

# Third party modules #
import sh

###############################################################################
class Assembly(object):
    """The co-assembly of all pools."""

    all_paths = """
    /output/
    """

    def __repr__(self): return '<%s object of %s>' % (self.__class__.__name__, self.parent)

    def __init__(self, aggregate):
        # Save parent #
        self.parent, self.aggregate = aggregate, aggregate
        # Auto paths #
        self.base_dir = self.parent.p.assembly_dir
        self.p = AutoPaths(self.base_dir, self.all_paths)

    def assemble(self):
        pairs = [(p.cleaner.fwd, p.cleaner.rev) for p in self.parent]
        sh.mpiexec('-n', 1, '-k', 81, '-o', self.p.output, *flatten([('-p', f, r) for f,r in pairs]))

    def make_plots(self):
        for cls_name in assembly_plots.__all__:
            cls = getattr(assembly_plots, cls_name)
            cls(self).plot()