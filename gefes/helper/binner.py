# Futures #
from __future__ import division

# Built-in modules #
from math import log10

# Internal modules #
from gefes.common.autopaths import AutoPaths
from gefes.helper.clusterer import Clusterer
from gefes.common.cache import property_cached
from gefes.fasta.single import FASTA
from gefes.helper.contig import Contig

# Third party modules #
from pandas import DataFrame

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
        self.bins_clusterer = None

    def run(self):
        if self.clusterer.clusters is None: self.clusterer.run()
        c_ids=list(set([c[1] for c in self.clusterer.clusters]))
        self.bins = dict.fromkeys(c_ids)
        for k in self.bins: self.bins[k] = []
        for co,cl in self.clusterer.clusters:
            self.bins[cl].append(co)
        self.bins_clusterer = self.clusterer.default


    def linkage_quality(self):
        linkage=[p.mapper.linkage for p in self.aggregate]
        
        qualities = {}
        all_contigs = list(self.bins_clusterer.frame.index)
        
        for k,b in self.bins.iteritems():
            qualities[k] = [0,0,0,0]
            for contig1 in b:
                for contig2 in b:
                    temp = 0
                    for p in linkage:
                        temp = temp + (sum(p[contig1][contig2]) > 0)
                        qualities[k][2] = qualities[k][2] + sum(p[contig1][contig2])
                    if temp > 0:
                        qualities[k][0] = qualities[k][0] + 1
                not_contigs = [c for c in all_contigs if c not in b]
                for contig2 in  not_contigs:
                    temp = 0
                    for p in linkage:
                        temp = temp + (sum(p[contig1][contig2]) > 0)
                        qualities[k][3] = qualities[k][3] + sum(p[contig1][contig2])
                    if temp > 0:
                        qualities[k][1] = qualities[k][1] + 1
                    
        return qualities
                
                    
    def bins_stats(self):
        out={}
        for k,b in self.bins.iteritems():
            out[k]=self.bins_clusterer.frame[[c for c in self.bins_clusterer.frame if "pool" in c]].loc[b].applymap(lambda x: log10(1+x))
            out[k]=out[k].median()
        return DataFrame.from_dict(out)
