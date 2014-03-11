# Futures #
from __future__ import division

# Built-in modules #
from itertools import product

# Internal modules #
from gefes.common.autopaths import AutoPaths
from gefes.common.cache import property_cached
from Bio.Seq import Seq

# Third party modules #

###############################################################################
class Contig(object):
    """Has a nucleotide frequency."""

    all_paths = """
    /lorem.txt
    /genes/ 
    """

    tetramers = ["".join(tetramer) for tetramer in product('ACGT', repeat=4)]
    for t in tetramers:
        if Seq(t).reverse_complement().tostring() != t:
            tetramers.remove(Seq(t).reverse_complement().tostring())
        for t in tetramers:
            if t[::-1] in tetramers:
                if t!=t[::-1]:
                    tetramers.remove(t[::-1])
    tetra_cats = {t:list(set((t,t[::-1],Seq(t).reverse_complement().tostring(),Seq(t[::-1]).reverse_complement().tostring()))) for t in tetramers}



    
    def __repr__(self): return '<%s object of %s>' % (self.__class__.__name__, self.parent)

    def __init__(self, assembly, record):
        # Save parent #
        self.parent, self.assembly = assembly, assembly
        # Auto paths #
        self.base_dir = self.parent.base_dir
        self.p = AutoPaths(self.base_dir, self.all_paths)
        self.record = record
        self.name = self.record.id
        
    def get_nuc_freq(self, windowsize):
        """Returns frequency of nucelotide sequence with length windowsize in
        given sequence."""
        freqs = {}

        # only A, C, G and T are allowed nucleotide characters
        allowed_nucs = {'A': 0, 'C': 0, 'G': 0, 'T': 0}
        is_nuc = lambda x: x in allowed_nucs

        upper = str(self.record.seq).upper()
        i = 0

        while i < len(upper) - windowsize + 1:
            nuc_win = upper[i:i + windowsize]

            # check if all nucleotides in the window are A,C,G or T
            for j, n in enumerate(nuc_win):
                if not is_nuc(n):
                    i += j + 1
                    break
            else:
                # if no break add nuc_win count
                freqs[nuc_win] = freqs.get(nuc_win, 0) + 1
                i += 1

        # normalize counts by number of windows to get frequency
        for nuc_win in freqs:
            freqs[nuc_win] = freqs[nuc_win] / float(self.length - windowsize + 1)
        return freqs

    def get_all_tetra_nuc_freqs(self):
        freqs=[]
        for t in Contig.tetra_cats:
            freqs = freqs + [sum([self.tetra_nuc_freq.get(tt,0) for tt in Contig.tetra_cats[t]])]
        return freqs

    @property
    def gc_content(self):
        outp = 0
        gs = self.nuc_freq.get("G")
        cs = self.nuc_freq.get("C")
        if gs:
            outp = outp + gs
        if cs:
            outp = outp + cs
        return outp
    
    @property
    def length(self):
        return len(self.record.seq)

    @property_cached
    def nuc_freq(self):
        return self.get_nuc_freq(1)

    @property_cached
    def tetra_nuc_freq(self):
        return self.get_nuc_freq(4)

