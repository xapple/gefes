# Futures #
from __future__ import division

# Built-in modules #

# Internal modules #
from gefes.common.autopaths import AutoPaths
from gefes.helper.assembler import Assembly
from gefes.helper.binner import Binner
from gefes.graphs import aggregate_plots

###############################################################################
class Aggregate(object):
    """A arbitrary aggregate of several pools."""

    all_paths = """
    /logs/
    /graphs/
    /assembly/
    /binning/
    """

    def __repr__(self): return '<%s object "%s" with %i pools>' % \
                               (self.__class__.__name__, self.name, len(self))
    def __iter__(self): return iter(self.pools)
    def __len__(self): return len(self.pools)
    def __getitem__(self, key):
        if isinstance(key, basestring): return [c for c in self.children if str(c) == key][0]
        elif isinstance(key, int): return self.children[key]
        else: raise TypeError('key')

    @property
    def first(self): return self.pools[0]

    def run_pools(self, steps=None, **kwargs):
        for p in self.pools: p.runner.run()

    def run_pools_slurm(self, steps=None, **kwargs):
        return [p.run_slurm(steps, **kwargs) for p in self.pools]

    def __init__(self, name, pools, out_dir):
        # Attributes #
        self.name = name
        self.pools = pools
        self.out_dir = out_dir
        # Dir #
        self.base_dir = self.out_dir + self.name + '/'
        self.p = AutoPaths(self.base_dir, self.all_paths)
        # Common init stuff #
        self.load()

    def load(self):
        # Children #
        self.assembly = Assembly(self)
        self.binner = Binner(self)
        # All the plots #
        self.graphs = [getattr(aggregate_plots, cls_name)(self) for cls_name in aggregate_plots.__all__]

    def assemble(self):
        self.assembly.assemble()

    def make_plots(self):
        for graph in self.graphs: graph.plot()