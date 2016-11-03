# Built-in modules #
from collections import OrderedDict

# Internal modules #
from gefes.assemble.contig import Contig
from gefes.binning.concoct import Concoct
from gefes.outputs.hit_profile       import HitProfile
from gefes.outputs.trait_annotations import TraitAnnotations
from gefes.outputs.bins_summary      import BinsSummary
from gefes.report.assembly           import AssemblyReport

# First party modules #
from fasta import FASTA
from plumbing.cache import property_cached
from plumbing.autopaths import AutoPaths

# Third party modules #
import pandas

###############################################################################
class Assembler(object):
    """Inherit from this."""

    def __nonzero__(self): return bool(self.p.filtered)
    def __repr__(self): return '<%s object kmer %i>' % (self.__class__.__name__, self.kmer_size)
    def __len__(self):  return len(self.samples)
    def __getitem__(self, key):
        if isinstance(key, basestring): return [c for c in self.children if c.name == key][0]
        return self.children[key]

###############################################################################
class AssemblyResults(object):
    """Inherit from this."""

    def __nonzero__(self): return bool(self.contigs_fasta)
    def __iter__(self): return iter(self.contigs)
    def __init__(self, parent):
        self.parent = parent
        self.p      = parent.p
        self.contigs_fasta = FASTA(self.parent.p.filtered)

    @property_cached
    def contigs(self):
        """All the contigs produced returned as a list of our Contig custom objects.
        TODO: contig numbers and ids have an often of 1"""
        return [Contig(self.parent, record, num=i) for i,record in enumerate(self.contigs_fasta)]

    @property_cached
    def contig_id_to_contig(self):
        """A dictionary with contig names as keys and contig objects as values."""
        return {c.name: c for c in self.contigs}

    @property_cached
    def total_bp(self):
        """The total amount of base pairs generated in this assembly."""
        return sum(self.contigs_fasta.lengths)

    @property_cached
    def mappings(self):
        """Map each of the samples used in the assembly back to this assembly.
        TODO: This should be updated to use a directory in the assembly results directory
        and to remove the attributes from the Sample objects."""
        return OrderedDict([(s.name, getattr(s, "mapper_%i" % self.parent.kmer_size)) for s in self.parent.samples])

    @property_cached
    def binner(self):
        """Put the contigs of this assembly into bins."""
        return Concoct(self.parent.samples, self.parent, self.parent.p.bins_dir)

    @property_cached
    def mappings_per_sample(self):
        """A dataframe with columns as samples and rows as contig names.
        Values are mapped read counts."""
        get_counts = lambda s: s.mappers[self.parent].results.raw_contig_counts
        frame = pandas.DataFrame({s.name: get_counts(s) for s in self.parent.samples})
        return frame

    #----------------------------------- Outputs -----------------------------#
    @property_cached
    def hit_profile(self):
        """Gather all the information to make a profile of the hits from a search."""
        return HitProfile(self.parent, self.p.hit_profile_dir)

    @property_cached
    def trait_annotations(self):
        """Special zooming into a particular pathway."""
        return TraitAnnotations(self.parent, self.p.traits_dir)

    @property_cached
    def bins_summary(self):
        """Special summary of bin properties."""
        return BinsSummary(self.parent, self.p.bins_summary_dir)
