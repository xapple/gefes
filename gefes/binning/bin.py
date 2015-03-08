# Built-in modules #

# Internal modules #
from gefes.annotation.prokka import Prokka
from gefes.annotation.cogs   import SingleCOGs
from gefes.annotation.checkm import Checkm

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
    /evaluation/
    """

    def __repr__(self): return '<%s object "%s">' % (self.__class__.__name__, self.name)

    def __init__(self, binner, contigs, result_dir=None, num=None, name=None):
        """You have to specify one of either a `name` or a `num`."""
        # Save Attributes #
        self.binner = binner
        self.contigs = contigs
        # Num #
        if num is None: self.num = 0
        else:           self.num = num
        # Name #
        if name is None: self.name = str(self.num)
        else:            self.name = name
        # Base directory #
        if result_dir is None: self.result_dir = self.binner.p.bins_dir
        else:                  self.result_dir = result_dir
        # Auto paths #
        self.base_dir = self.result_dir + self.name + '/'
        self.p = AutoPaths(self.base_dir, self.all_paths)

    @property_cached
    def fasta(self):
        """A fasta file containing only the contigs pertaining to this bin in it."""
        fasta = FASTA(self.p.fasta)
        if not fasta.exists:
            with fasta as handle:
                for contig in self.contigs:
                    handle.add_seq(contig.record)
        return fasta

    @property_cached
    def annotation(self):
        """The results from annotation the bin."""
        return Prokka(self, self.p.annotation_dir)

    @property_cached
    def single_cogs(self):
        """The results from finding single copy COGs in the bin."""
        return SingleCOGs(self, self.p.annotation_dir)

    @property_cached
    def evaluation(self):
        """The results from evaluating the bin completness."""
        return Checkm(self, self.p.evaluation_dir)