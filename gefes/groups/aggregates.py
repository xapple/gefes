# Futures #
from __future__ import division

# Built-in modules #
import gefes
from gefes.assemble.ray             import Ray
from gefes.binning.concoct          import Concoct
from gefes.running.aggregate_runner import AggregateRunner
from gefes.report.aggregate         import AggregateReport

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
        if isinstance(key, basestring): return [c for c in self.children if c.name == key][0]
        elif isinstance(key, int): return self.children[key]
        elif isinstance(key, slice): return self.children[key]
        else: raise TypeError('key')

    def __init__(self, name, samples, base_dir=None):
        # Attributes #
        self.name = name
        self.samples, self.children = samples, samples
        # Base directory #
        if base_dir == None: self.base_dir = gefes.view_dir + 'aggregates/' + name + '/'
        else: self.base_dir = base_dir
        # Delayed init #
        self.loaded = False

    def load(self):
        """A delayed kind of __init__ that is not called right away to avoid
        crowding the RAM of the python interpreter when you just import gefes"""
        # Load #
        self.loaded = True
        # Automatic paths #
        self.p = AutoPaths(self.base_dir, self.all_paths)
        # Assemble #
        self.assembly = Ray(self.samples, self.p.assembly_dir)
        # Binner #
        self.binner = Concoct(self.samples, self.assembly, self.p.binning_dir)
        # Annotation #
        #self.phylotyper = Phylotyper(self)
        #self.annotation = Prokka(self)
        # Runner #
        self.runner = AggregateRunner(self)
        # Report #
        self.report = AggregateReport(self)
        # For convenience #
        return self

    @property
    def first(self): return self.samples[0]