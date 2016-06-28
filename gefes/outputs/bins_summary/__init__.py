# Built-in modules #

# Internal modules #
import os, inspect
from collections import defaultdict

# First party modules #
from plumbing.cache     import property_cached
from plumbing.autopaths import AutoPaths
from fasta import FASTA
from seqsearch import SeqSearch
from seqsearch.blast import BLASTdb

# Third party modules #
import pandas

# Current directory #
filename    = inspect.getframeinfo(inspect.currentframe()).filename
current_dir = os.path.dirname(os.path.abspath(filename)) + '/'

###############################################################################
class BinsSummary(object):
    """Let's make some summaries."""

    short_name = "bins_summary"

    all_paths = """
    /bins_to_gc.tsv
    """

    def __nonzero__(self): return bool(self.p.hits)

    def __init__(self, assembly, base_dir):
        # Save parent #
        self.parent, self.assembly = assembly, assembly
        self.binner                = self.assembly.results.binner
        self.contig_to_bin         = self.binner.results.contig_id_to_bin_id
        # Auto paths #
        self.base_dir = base_dir
        self.p = AutoPaths(self.base_dir, self.all_paths)

    def run(self, verbose=True):
        # Make the TSVs #
        self.bins_to_gc.to_csv(self.p.bins_to_gc.path, sep='\t', float_format='%.5g')
        # Return #
        return self.p.traits_x_bins

    @property_cached
    def bins_to_gc(self):
        """A matrix with linking every bin to it's total average GC content.
        Note: you can't just average the values of the contigs average GC,
        it has to be weighted by length. Or just contact all contigs, and
        calculate GC on that."""
        # Build a new frame #
        result = defaultdict(float)
        for bin in self.binner.results.bins:
            result[bin.name] = bin.average_gc
        # Return #
        result = pandas.Series(result)
        return result

    @property_cached
    def results(self):
        results = BinsSummaryResults(self)
        message = "You can't access results from the bins summary before running the algorithm."
        if not results: raise Exception(message)
        return results

###############################################################################
class BinsSummaryResults(object):

    def __nonzero__(self): return bool(self.p.hits)

    def __init__(self, bs):
        self.bs, self.parent = bs, bs
        self.p = bs.p

    @property_cached
    def bins_to_gc(self):
        return pandas.io.parsers.read_csv(self.p.bins_to_gc.path,
                                          sep='\t', index_col=0, encoding='utf-8')