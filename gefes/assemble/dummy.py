# Futures #
from __future__ import division

# Built-in modules #

# Internal modules #
from gefes.assemble import Assembler, AssemblyResults

# First party modules #
from plumbing.autopaths import AutoPaths
from plumbing.cache import property_cached

# Third party modules #

###############################################################################
class DummyAssembler(Assembler):
    """For when it's already been done by someone else."""

    short_name = 'megahit'
    long_name  = 'MEGAHIT assembler v1.0.6.1'
    executable = 'megahit'
    url        = 'https://github.com/voutcn/megahit'
    dependencies = []

    all_paths = Assembler.all_paths + """
    /output/final.contigs.fa
    /output/log
    /output/opts.txt
    /filtered.fasta
    /stdout.txt
    /stderr.txt
    """

    kmer_size = 'variable'

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
        return "Megahit with kmer %s" % self.kmer_size

    @property
    def description(self):
        return "Megahit with kmer %s and %s bp cutoff" % (self.kmer_size, self.length_cutoff)

    @property_cached
    def results(self):
        results = MegahitResults(self)
        if not results: raise Exception("You can't access results from Megahit before running the assembly.")
        return results

###############################################################################
class DummyAssemblerResults(AssemblyResults):
    pass