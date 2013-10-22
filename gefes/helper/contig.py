# Futures #
from __future__ import division

# Built-in modules #

# Internal modules #
from gefes.common.autopaths import AutoPaths

# Third party modules #

###############################################################################
class Contig(object):
    """Has a nucleotide frequency, a coverage, etc.."""

    all_paths = """
    /lorem.txt
    """

    def __repr__(self): return '<%s object of %s>' % (self.__class__.__name__, self.parent)

    def __init__(self, assembly, record):
        # Save parent #
        self.parent, self.assembly = assembly, assembly
        # Auto paths #
        self.base_dir = self.parent.base_dir
        self.p = AutoPaths(self.base_dir, self.all_paths)
        self.record = record

    def get_nuc_freq(sequence):
        """Returns character frequency of given string sequence."""
        freqs = {}

        upper = sequence.upper()
        for c in upper:
            freqs[c] = freqs.get(c, default=0) + 1
