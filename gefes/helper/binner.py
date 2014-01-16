# Futures #
from __future__ import division

# Built-in modules #
from itertools import product

# Internal modules #
from gefes.common.autopaths import AutoPaths
from gefes.helper.clusterer import Clusterer
from gefes.common.cache import property_cached
from Bio.Seq import Seq
from gefes.fasta.single import FASTA


# Third party modules #
import pandas

###############################################################################
class Binner(object):
    """Makes the matrix detailing all information for every contig,
    The output is N bins, each with one or more contigs inside."""

    all_paths = """
    /clustering/
    /bins/
    /frame.csv
    """

    def __repr__(self): return '<%s object of %s>' % (self.__class__.__name__, self.parent)

    def __init__(self, aggregate):
        # Save parent #
        self.parent, self.aggregate = aggregate, aggregate
        # Auto paths #
        self.base_dir = self.parent.p.binning_dir
        self.p = AutoPaths(self.base_dir, self.all_paths)
        # Children #
        self.clusterer=Clusterer(self)
        # Output #
        self.bins = []

    @property_cached
    def frame(self):
        tetramers = ["".join(tetramer) for tetramer in product('ACGT', repeat=4)]
        for t in tetramers:
            if Seq(t).reverse_complement().tostring() != t:
                tetramers.remove(Seq(t).reverse_complement().tostring())
        for t in tetramers:
            if t[::-1] in tetramers:
                if t!=t[::-1]:
                    tetramers.remove(t[::-1])
        tetra_cats={t:list(set((t,t[::-1],Seq(t).reverse_complement().tostring(),Seq(t[::-1]).reverse_complement().tostring()))) for t in tetramers}
        columns = ['length'] + [s.id_name for s in self.aggregate] + ['freq_' + t for t in tetramers]
        rows = [c.name for c in self.aggregate.assembly.contigs]
        data = [[c.length] +
                [s.mapper.coverage[c.name]["cov_mean"] for s in self.aggregate] +
                [sum([c.tetra_nuc_freq.get(tt,0) for tt in tetra_cats[t]]) for t in tetra_cats] for c in self.aggregate.assembly.contigs]
        return pandas.DataFrame(data, columns=columns, index=rows)

    def filtered_frame(self,max_freq,min_len):
        temp_frame = self.frame
        if(min_len is not None):
            temp_frame = temp_frame[temp_frame.length > min_len]
        if(max_freq is not None):
            good_ones = temp_frame[[c for c in temp_frame if "freq" in c]].apply(lambda x: sum(x>max_freq)==0,1)
            temp_frame = temp_frame[good_ones]
        return temp_frame

    
    def export_frame(self):
        self.frame.to_csv(self.p.frame)

