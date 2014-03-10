# Futures #
from __future__ import division

# Built-in modules #
from math import log10
import os

# Internal modules #
from gefes.common.autopaths import AutoPaths
from gefes.helper.clusterer import Clusterer
from gefes.common.cache import property_cached
from gefes.fasta.single import FASTA
from gefes.helper.contig import Contig
from gefes.common.autopaths import AutoPaths

# Third party modules #
from pandas import DataFrame

###############################################################################
class Binner(object):
    """Makes the matrix detailing all information for every contig,
    The output is N bins, each with one or more contigs inside."""

    all_paths = """
    """

    def __repr__(self): return '<%s object of %s with %i Binnings>' % \
                               (self.__class__.__name__, self.parent, len(self))
    def __iter__(self): return iter(self.binnings)
    def __len__(self): return len(self.binnings)
    def __getitem__(self, key):
        return self.binnings[key]
            
    def __init__(self, parent):
        # Save parent #
        self.parent = parent
        self.assembly = parent.assembly
        # Auto paths #
        self.base_dir = self.parent.p.binning_dir
        self.p = AutoPaths(self.base_dir, self.all_paths)
        # Children #
        self.binnings=self.load()


    def load(self):
        binnings={}
        binnings_paths = os.listdir(self.base_dir)
        for b in binnings_paths:
            print b
            binnings[b]=Binning(self, b)
        return binnings

    def new(self,name,clusterer):
        self.binnings[name] = Binning(self, name, clusterer)
             

###############################################################################        
class Binning(object):

    all_paths = """
    /bins/
    /settings.json
    /linkage_quality.csv
    /linkage_matrix.csv
    /coverage_stats.csv
    /frame.csv
    """

    def __init__(self, parent,name ,clusterer = None):
        # Save parent #
        self.parent =  parent
        self.name = name
        # Auto paths #
        self.base_dir = self.parent.p._base_dir+name
        self.p = AutoPaths(self.base_dir, self.all_paths)
        # Output #
        self.clusterer = clusterer
        if clusterer is None: self.load()

    def load(self):
        pass
    
            
    def run(self):
        if self.clusterer is None: raise Exception('Clustering has already been run, or is broken, or not existing')
        self.clusterer.run(self.parent.assembly)
        c_ids=list(set([c[1] for c in self.clusterer.clusters]))
        self.bins = dict.fromkeys(c_ids)
        for k in self.bins: self.bins[k] = []
        for co,cl in self.clusterer.clusters:
            self.bins[cl].append(co)
        self.frame=self.clusterer.frame
        self.export()
            
    def linkage_quality(self):
        linkage=[p.mapper.linkage for p in self.aggregate]
        
        qualities = {}
        all_contigs = list(self.frame.index)
        
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

    def linkage_matrix(self):
        linkage=[p.mapper.linkage for p in self.aggregate]
        
        matrix = DataFrame(0,index=self.bins.keys(),columns=self.bins.keys())
        for k1,b1 in self.bins.iteritems():
            for k2,b2 in self.bins.iteritems():
                    for contig1 in b1:
                        for contig2 in b2:
                            temp = 0
                            for p in linkage:
                                if contig1 != contig2:
                                    matrix[k1][k2] = matrix[k1][k2] + sum(p[contig1][contig2])/2.0
                                    matrix[k2][k1] = matrix[k2][k1] + sum(p[contig1][contig2])/2.0
        return matrix
                    
    def bins_stats(self):
        out={}
        for k,b in self.bins.iteritems():
            out[k]=self.frame[[c for c in self.frame if "pool" in c]].loc[b].applymap(lambda x: log10(1+x))
            out[k]=out[k].median()
        return DataFrame.from_dict(out)


    def export(self):
        path=self.p.bins
        for name,contig_list in self.bins.iteritems():
            self.parent.assembly.write_contiglist(contig_list,path,"bin_"+str(name)+".fasta")

