# Futures #
from __future__ import division

# Built-in modules #

# Internal modules #
from plumbing.autopaths import AutoPaths

###############################################################################
class Aggregate(object):
    """A arbitrary aggregate of several samples."""

    all_paths = """
    /logs/
    /graphs/
    /assembly/
    /binning/
    """

    def __repr__(self): return '<%s object "%s" with %i samples>' % \
                               (self.__class__.__name__, self.name, len(self))
    def __iter__(self): return iter(self.samples)
    def __len__(self): return len(self.samples)
    def __getitem__(self, key):
        if isinstance(key, basestring): return [c for c in self.children if str(c) == key][0]
        elif isinstance(key, int): return self.children[key]
        else: raise TypeError('key')

    @property
    def first(self): return self.samples[0]

    def run_samples(self, steps=None, **kwargs):
        for p in self.samples: p.runner.run()

    def run_samples_slurm(self, steps=None, **kwargs):
        return [p.run_slurm(steps, **kwargs) for p in self.samples]

    def __init__(self, name, samples, out_dir):
        # Attributes #
        self.name = name
        self.samples = samples
        self.out_dir = out_dir
        # Dir #
        self.base_dir = self.out_dir + self.name + '/'
        self.p = AutoPaths(self.base_dir, self.all_paths)

    def assemble(self):
        self.assembly.assemble()

    def make_plots(self):
        for graph in self.graphs: graph.plot()