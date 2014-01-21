# Futures #
from __future__ import division

# Built-in modules #

# Internal modules #
from gefes.common.autopaths import AutoPaths
from gefes.common.cache import property_cached

# Third party modules #
import os
import pandas
from sklearn import cluster
from scipy.spatial import distance
from scipy.cluster import hierarchy
from sklearn.cluster  import KMeans
import math

###############################################################################
class Clusterer(object):
    """Recieves the matrix detailing all information for every contig,
    and is responsible for deciding which contigs go together."""

    all_paths = """
    /lorem.txt
    /kmeans/
    """

    def __repr__(self): return '<%s object of %s>' % (self.__class__.__name__, self.parent)

    def __init__(self, binner):
        # Save parent #
        self.parent = binner
        self.assembly = self.parent.aggregate.assembly
        # Auto paths #
        self.base_dir = self.parent.p.clustering
        self.p = AutoPaths(self.base_dir, self.all_paths)
        # Methods
        self.kmeans = GefesKMeans(self)
 
    def run(self):
        self.kmeans.run()
        
 
###############################################################################
class GefesKMeans(object):
    """Receives the matrix and uses kmeans to cluster it in N different (linearly separated) clusters"""

    all_paths = """
    /tetramer_clusters/
    /coverage_clusters/
    """

    
    def __init__(self,parent,nb = 8,max_freq = None,min_length = None):
        # Save parent
        self.parent = parent
        # Kmeans and params
        self.number_clusts=nb
        self.algorithm = KMeans(self.number_clusts)
        # Autopath
        self.base_dir = self.parent.p.kmeans
        self.p = AutoPaths(self.base_dir, self.all_paths)
        # Clusters
        self.tetras_clusters = None
        self.coverage_clusters = None
        # Filters
        self.max_freq = max_freq
        self.min_length = min_length
        self.log10covs=False
        

    def run(self, nb = None, max_freq = None, min_len = None,log10covs=None):
        if nb:
            self.number_clusts=nb
            self.algorithm = KMeans(self.number_clusts)
        if max_freq: self.max_freq = max_freq
        if min_len: self.min_length = min_len
        if log10covs: self.log10covs = log10covs
        self.frame = self.parent.assembly.filtered_frame(self.max_freq,self.min_length)
        tetras = self.frame[[c for c in self.frame if "freq" in c]]
        covers = self.frame[[c for c in self.frame if "pool" in c]]
        log10p1 = lambda x: math.log10(x+1)
        if self.log10covs: covers = covers.applymap(log10p1)
        self.tetras_clusters = self.algorithm.fit_predict(tetras)        
        self.coverage_clusters = self.algorithm.fit_predict(covers)
        self.export()

    def export(self):
        parameters = "-".join(["nb_clusts",str(self.number_clusts),"min_len",str(self.min_length),"max_freq",str(self.max_freq)])
        assembly = self.parent.parent.parent.assembly
        contigs = [f[0] for f in self.frame.itertuples()]
        path=os.path.join(self.p.tetramer,parameters)
        for i in set(self.tetras_clusters):
            contig_list=[]
            for j in range(0,len(contigs)):
                if(self.tetras_clusters[j] == i): contig_list.append(contigs[j])
            assembly.write_contiglist(contig_list,path,"bin_"+str(i)+".fasta")
        path=os.path.join(self.p.coverage,parameters)
        for i in set(self.coverage_clusters):
            contig_list=[]
            for j in range(0,len(contigs)):
                if(self.coverage_clusters[j] == i): contig_list.append(contigs[j])
            assembly.write_contiglist(contig_list,path,"bin_"+str(i)+".fasta")

            
