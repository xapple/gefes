# Built-in modules #

# Internal modules #
from gefes.annotation.prokka import Prokka

# First party modules #
from plumbing.autopaths import AutoPaths
from plumbing.cache import property_cached
from fasta import FASTA

# Third party modules #

###############################################################################
class Bin(object):
    """A Bin is a collection of Contigs that were identified as potentially
    coming from the same population/species/strain of organisms."""

    all_paths = """
    /contigs.fasta
    /annotation/
    """

    def __repr__(self): return '<%s object "%s">' % (self.__class__.__name__, self.name)

    def __init__(self, binner, contigs, result_dir=None, num=None, name=None):
        # Save Attributes #
        self.binner = binner
        self.contigs = contigs
        self.num = int(num)
        # Base directory #
        if result_dir is None: self.result_dir = self.binner.p.bins_dir
        else:                  self.result_dir = result_dir
        # Name #
        if name is None: self.name = str(self.num)
        else:            self.name = name
        # Auto paths #
        self.base_dir = self.result_dir + self.name + '/'
        self.p = AutoPaths(self.base_dir, self.all_paths)
        # Extra objects #
        self.annotation = Prokka(self, self.p.annotation_dir)
        #self.reassembly = Ray()

    @property_cached
    def fasta(self):
        """A fasta file containing only the contigs pertaining to this bin in it."""
        fasta = FASTA(self.p.fasta)
        if not fasta.exists:
            with fasta as handle:
                for contig in self.contigs:
                    handle.add_seq(contig.record)
        return fasta
