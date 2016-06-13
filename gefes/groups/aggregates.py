# Futures #
from __future__ import division

# Built-in modules #
from collections import OrderedDict

# Internal modules #
import gefes
from gefes.assemble.ray             import Ray
from gefes.running.aggregate_runner import AggregateRunner
from gefes.report.aggregate         import AggregateReport
from gefes.merged.newbler           import Newbler
from gefes.status.projects          import ProjectStatus

# First party modules #
from plumbing.autopaths import AutoPaths

###############################################################################
class Aggregate(object):
    """A arbitrary aggregate of several samples.
    Typically, a `Project` object will inherit from this
    and extent the load() method."""

    all_paths = """
    /logs/
    /graphs/
    /assembly/
    /merged/
    /report/report.pdf
    """

    def __repr__(self): return '<%s object "%s" with %i samples>' % \
                               (self.__class__.__name__, self.name, len(self))
    def __iter__(self): return iter(self.samples)
    def __len__(self): return len(self.samples)
    def __getitem__(self, key):
        if isinstance(key, basestring): return [c for c in self.children if c.name == key][0]
        return self.children[key]

    def __init__(self, name, samples, base_dir=None):
        # Attributes #
        self.name = name
        self.samples, self.children = samples, samples
        # Base directory #
        if base_dir is None: self.base_dir = gefes.view_dir + 'aggregates/' + name + '/'
        else: self.base_dir = base_dir
        # Delayed init #
        self.loaded = False

    def load(self):
        """A delayed kind of __init__ that is not called right away to avoid
        crowding the RAM of the python interpreter when you just import gefes."""
        # Load #
        self.loaded = True
        # Automatic paths #
        self.p = AutoPaths(self.base_dir, self.all_paths)
        # Assemble #
        self.assembly = Ray(self.samples, self.p.assembly_dir)
        # Assemble with different kmer sizes #
        self.assembly_51 = Ray(self.samples, self.p.assembly_dir, kmer_size=51)
        self.assembly_61 = Ray(self.samples, self.p.assembly_dir, kmer_size=61)
        self.assembly_71 = self.assembly
        self.assembly_81 = Ray(self.samples, self.p.assembly_dir, kmer_size=81)
        # Also combine the different kmer sizes #
        self.merged = Newbler(self.samples, self.assemblies.values(), self.p.merged_dir)
        # Annotation #
        #self.phylotyper = Phylotyper(self)
        #self.annotation = Prokka(self)
        # Runner #
        self.runner = AggregateRunner(self)
        # Report #
        self.report = AggregateReport(self)
        # Status #
        self.status = ProjectStatus(self)
        # For convenience #
        return self

    @property
    def first(self): return self.samples[0]

    @property
    def assemblies(self):
        """A dictionary useful for trying different assemblies of different sizes.
        Keys are kmer-sizes and values are assembler objects"""
        return OrderedDict(((51, self.assembly_51),
                            (61, self.assembly_61),
                            (71, self.assembly_71),
                            (81, self.assembly_81)))