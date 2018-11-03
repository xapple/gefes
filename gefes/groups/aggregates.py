# Futures #
from __future__ import division

# Built-in modules #
import warnings
from collections import OrderedDict

# Internal modules #
from gefes.assemble.ray             import Ray
from gefes.assemble.megahit         import Megahit
from gefes.assemble.dummy           import DummyAssembler
from gefes.running.aggregate_runner import AggregateRunner
from gefes.report.aggregate         import AggregateReport
from gefes.merged.newbler           import Newbler
from gefes.status.projects          import ProjectStatus

# First party modules #
from plumbing.autopaths import AutoPaths, DirectoryPath
from plumbing.cache     import property_cached

###############################################################################
class Aggregate(object):
    """A arbitrary aggregate of several samples.
    Typically, a `Project` object will inherit from this
    and extent the load() method."""

    default_assembler = "megahit"

    all_paths = """
    /samples/
    /logs/
    /graphs/
    /assembly/
    /merged/
    /report/report.pdf
    """

    def __repr__(self): return '<%s object "%s" with %i samples>' % \
                               (self.__class__.__name__, self.name, len(self))

    def __str__(self):             return self.short_name
    def __iter__(self):            return iter(self.children)
    def __len__(self):             return len(self.children)
    def __contains__(self, item):  return item in self.children

    def __add__(self, other):
        name     = self.name + " and " + other.name
        children = self.children + other.children
        return self.__class__(name, children, sort=False)

    def __getitem__(self, key):
        if   isinstance(key, basestring):  return [c for c in self.children if str(c) == key][0]
        elif isinstance(key, int):
            if key < 0:                    return self.children[key]
            if hasattr(self.first, 'num'): return [c for c in self.children if int(c.num) == key][0]
            else:                          return self.children[key]
        elif isinstance(key, slice):       return self.children[key]
        else:                              raise TypeError('key')

    @property
    def first(self):  return self.children[0]
    @property
    def second(self): return self.children[1]
    @property
    def third(self):  return self.children[2]

    def __init__(self, name, samples, out_dir, sort=True):
        # Attributes #
        self.name       = name
        self.short_name = name
        self.samples    = samples
        self.children   = samples
        # Check names are unique #
        names = [s.short_name for s in self.samples]
        assert len(names) == len(set(names))
        # Are the samples numbered #
        have_numbers = all(s.info.get('sample_num') for s in samples)
        if not have_numbers: warnings.warn("Not all samples of project '%s' were numbered." % self)
        # Sort the samples #
        if sort:
            if have_numbers: samples.sort(key=lambda s: int(s.info['sample_num']))
            else:            samples.sort(key=lambda s: s.short_name)
        # Base directory #
        self.base_dir = DirectoryPath(out_dir + self.name + '/')

    #-------------------------------- Properties -----------------------------#
    @property_cached
    def p(self): return AutoPaths(self.base_dir, self.all_paths)

    @property_cached
    def assembly_51(self): return Ray(self.samples, self.p.assembly_dir, kmer_size=51)
    @property_cached
    def assembly_61(self): return Ray(self.samples, self.p.assembly_dir, kmer_size=61)
    @property_cached
    def assembly_71(self): return Ray(self.samples, self.p.assembly_dir, kmer_size=71)
    @property_cached
    def assembly_81(self): return Ray(self.samples, self.p.assembly_dir, kmer_size=81)

    @property_cached
    def merged(self):
        """All assemblies merged into a bigger one."""
        choices = {'megahit': (Megahit,        (self.samples, self.p.merged_dir)),
                   'newbler': (Newbler,        (self.samples, self.assemblies.values(), self.p.merged_dir)),
                   'dummy':   (DummyAssembler, (self.samples, self.p.merged_dir))}
        cls, params = choices.get(self.default_assembler)
        return cls(*params)

    @property_cached
    def runner(self):
        """The runner object."""
        return AggregateRunner(self)

    @property_cached
    def report(self):
        """The PDF report."""
        return AggregateReport(self)

    @property_cached
    def status(self):
        """The status."""
        return ProjectStatus(self)

    #-------------------------------- Shortcuts -----------------------------#
    @property
    def assembly(self):
        """Convenience shortcut. By default the 71 kmer assembly."""
        return self.assembly_71

    @property_cached
    def assemblies(self):
        """A dictionary useful for trying different assemblies of different sizes.
        Keys are kmer-sizes and values are assembler objects"""
        return OrderedDict(((51, self.assembly_51),
                            (61, self.assembly_61),
                            (71, self.assembly_71),
                            (81, self.assembly_81)))