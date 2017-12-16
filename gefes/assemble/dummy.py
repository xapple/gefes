# Futures #
from __future__ import division

# Built-in modules #
import os

# Internal modules #
from gefes.assemble import Assembler, AssemblyResults

# First party modules #
from plumbing.autopaths import AutoPaths
from plumbing.cache import property_cached
from fasta import FASTA

# Third party modules #

###############################################################################
class DummyAssembler(Assembler):
    """For when it's already been done by someone else."""

    short_name = 'metaspades'
    long_name  = 'meta-SPAdes assembler v3.10.0'
    executable = 'metaspades.py'
    url        = 'https://github.com/ablab/spades'
    dependencies = []

    all_paths = Assembler.all_paths + """
    /contigs.fasta
    /filtered.fasta
    /params.txt
    /spades.log
    """

    kmer_size     = 'variable'
    length_cutoff = 2500

    def __repr__(self): return '<%s object on samples "%s">' % \
        (self.__class__.__name__, ','.join(map(lambda s: s.short_name, self.samples)))

    def __init__(self, samples, result_dir, length_cutoff=1000):
        # Base parameters #
        self.samples       = samples
        self.children      = samples
        self.result_dir    = result_dir
        self.length_cutoff = length_cutoff
        # Auto paths #
        self.base_dir = self.result_dir + self.short_name + '/'
        self.p = AutoPaths(self.base_dir, self.all_paths)

    @property
    def short_description(self):
        return "SPAdes version: 3.10.0"

    @property
    def description(self):
        return "SPAdes version: 3.10.0 and %i bp cutoff" % self.length_cutoff

    def run(self):
        # Check there is something #
        contigs = FASTA(self.p.contigs)
        if len(contigs) == 0: raise Exception("DummyAssembler found exactly 0 contigs in your dataset.")
        # Filter short contigs #
        filtered = FASTA(self.p.filtered)
        contigs.extract_length(new_path=filtered, lower_bound=self.length_cutoff)
        # Make indexes (used later) #
        if filtered and not os.path.exists(filtered + '.fai'): filtered.index_samtools()

    @property_cached
    def results(self):
        results = DummyAssemblerResults(self)
        return results

###############################################################################
class DummyAssemblerResults(AssemblyResults):
    pass