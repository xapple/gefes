# Futures #

# Built-in modules #

# Internal modules #
import gefes
from gefes.taxonomy.phylophlan import Phylophlan

# First party modules #
from plumbing.autopaths import AutoPaths
from plumbing.cache     import property_cached

###############################################################################
class Lump(object):
    """Put several Aggregates objects (or Projects objects) together."""

    all_paths = """
    /taxonomy/
    """

    def __repr__(self): return '<%s object "%s" with %i children>' % \
                               (self.__class__.__name__, self.name, len(self))
    def __iter__(self): return iter(self.children)
    def __len__(self): return len(self.children)
    def __getitem__(self, key):
        if isinstance(key, basestring): return [c for c in self.children if c.name == key][0]
        return self.children[key]

    def __init__(self, name, aggregates, base_dir=None):
        # Attributes #
        self.name = name
        self.aggregates, self.children = aggregates, aggregates
        # Base directory #
        if base_dir is None: self.base_dir = gefes.view_dir + 'lumps/' + name + '/'
        else: self.base_dir = base_dir
        # Automatic paths #
        self.p = AutoPaths(self.base_dir, self.all_paths)

    @property
    def first(self): return self.children[0]

    @property
    def samples(self):
        return [s for a in self.aggregates for s in a]

    #-------------------------------------------------------------------------#
    def run_phylophlan(self, cpus=None):
        """Special method to run a combined PhyloPhlAn on all samples
        coming from multiple projects. We are going to cheat a lot below."""
        self.phylophlan.run(cpus=cpus)

    @property_cached
    def phylophlan(self): return Phylophlan(self, self.p.taxonomy_dir, self.name)

    @property_cached
    def good_bins(self):
        # Patch the names #
        for proj in self.aggregates:
            for b in proj.merged.results.binner.results.good_bins:
                b.name = proj.name[-2:] + '_' + b.name
        # Return #
        return [b for proj in self.aggregates for b in proj.merged.results.binner.results.good_bins]

    @property_cached
    def results(self):
        return type('FakeResults', (object,), {'good_bins': self.good_bins})
