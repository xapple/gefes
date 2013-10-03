# Futures #
from __future__ import division

# Built-in modules #

# Internal modules #
from gefes.common import AutoPaths

###############################################################################
class Collection(object):
    """A collection of aggregates."""

    def __repr__(self): return 'Collection: %s' % (self.children)
    def __iter__(self): return iter(self.children)
    def __len__(self): return len(self.children)

    def __init__(self, children):
        self.children = children

    @property
    def first(self): return self.children[0]

    def __getitem__(self, key):
        if isinstance(key, basestring):
            return [c for c in self.children if c.name == key.lower()][0]
        elif isinstance(key, int):
            if hasattr(self.first, 'num'): return [c for c in self.children if c.num == key][0]
            else: return self.children[key]
        else: raise TypeError('key')

###############################################################################
class Aggregate(object):
    """A arbitrary aggregate of several pools."""

    all_paths = """
    /graphs/
    /logs/
    /results/slurm_report.csv
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

    def __init__(self, name, pools, out_dir):
        # Attributes #
        self.name = name
        self.pools = pools
        self.loaded = False
        # Dir #
        self.base_dir = out_dir + self.name + '/'
        self.p = AutoPaths(self.base_dir, self.all_paths)

    def run_pools(self, steps=None, **kwargs):
        for p in self.pools: p()

    def run_pools_slurm(self, steps=None, **kwargs):
        return [p.run_slurm(steps, **kwargs) for p in self.pools]