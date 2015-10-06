# Futures #
from __future__ import division

# Built-in modules #

# Internal modules #
from gefes.annotation.prokka import Prokka
from gefes.annotation.prodigal import Prodigal
from gefes.taxonomy.phylosift import Phylosift

# First party modules #
from fasta import FASTA
from plumbing.autopaths import AutoPaths
from plumbing.cache import property_cached

###############################################################################
class Contig(object):
    """A contig as predicted by the assembler. It has for instance a nucleotide frequency
    and annotations. Possibilities for analyzing a contig's content:
    * Prokka
         (protein calling and function assignment from contigs)
    * PhyloPhlAn
        (Input the the (unannotated) protein sequences of the input genomes, get taxonomic predictions,
        more for complete genomes or bins)
    * MetaPhlAn
        (Input the raw reads and obtain some visualization)
    * PhyloSift
        (Input short or long sequence and have them compared against marker sets
        of 40 genes from the tree domains, uses hmmalign and infernal)
    * Blobology
        (Similar to the GC-vs-Coverage plots we make)
    """

    all_paths = """
    /contig.fasta
    /annotation/
    /taxonomy/
    """

    def __repr__(self): return '<%s object "%s">' % (self.__class__.__name__, self.name)
    def __len__(self): return len(self.record)

    def __init__(self, assembly, record, num=None):
        # Save parent #
        self.parent, self.assembly = assembly, assembly
        # Basic attributes #
        self.record = record
        self.name = self.record.id
        self.num = int(num)
        # Auto paths #
        self.base_dir = self.parent.base_dir + "contigs/" + str(self.num) + '/'
        self.p = AutoPaths(self.base_dir, self.all_paths)

    @property
    def length(self):
        return len(self.record.seq)

    @property_cached
    def fasta(self):
        """A fasta file containing only this contig."""
        fasta = FASTA(self.p.contig_fasta)
        if not fasta.exists:
            fasta.create()
            fasta.add_seq(self.record)
            fasta.close()
        return fasta

    @property_cached
    def annotation(self):
        """The annotation that can be made on this contig."""
        return Prokka(self, self.p.annotation_dir)

    @property_cached
    def proteins(self):
        """The predicted proteins that this contig contains."""
        return Prodigal(self, self.p.annotation_dir)

    @property_cached
    def ribosomal_proteins(self):
        """The predicted ribosomal proteins this contig contains."""
        pass

    @property_cached
    def taxonomy(self):
        """The predicted taxonomic information associated with this contig."""
        return Phylosift(self, self.p.taxonomy_dir)

    #-------------------------------------------------------------------------#
    @property
    def gc_content(self):
        pass

    #-------------------------------------------------------------------------#
    def get_nuc_freq(self, windowsize=4):
        """Returns frequency of nucleotides in this contig with length windowsize"""
        freqs = {}
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

    @property_cached
    def tetra_nuc_freq(self):
        return self.get_nuc_freq(4)

