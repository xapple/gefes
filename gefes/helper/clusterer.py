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
import scipy.stats

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
        self.default = self.kmeans
        self.clusters = None
        
    def run(self):
        self.kmeans.run()
        self.clusters=kmeans.clusters
        
    def log10(x):
        return covers.applymap(math.log10)

    def log10p1(x):
        lambda x: math.log10(x+1)
        return covers.applymap(log10p1)

    def rank(x):
        return x.apply(scipy.stats.rankdata,0)

    transforms = {'log10' : log10, 'log10p1' : log10p1 , 'rank' : rank}
    
        
###############################################################################
class GefesKMeans(object):
    """Receives the matrix and uses kmeans to cluster it in N different (linearly separated) clusters"""

    all_paths = """
    /tetramer_clusters/
    /coverage_clusters/
    /combo_clusters/
    """

    
    def __init__(self,parent,nb = 8,max_freq = None,min_length = None, transform = None):
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
        self.combo_clusters = None
        # Filters
        self.max_freq = max_freq
        self.min_length = min_length
        self.transform = None
        self.clusters = None
        

    def run(self, nb = None, max_freq = None, min_len = None, transform = None):
        if nb:
            self.number_clusts=nb
            self.algorithm = KMeans(self.number_clusts)
        if max_freq: self.max_freq = max_freq
        if min_len: self.min_length = min_len
        if transform: self.transform = transform
        self.frame = self.parent.assembly.filtered_frame(self.max_freq,self.min_length)
        tetras = self.frame[[c for c in self.frame if "freq" in c]]
        covers = self.frame[[c for c in self.frame if "pool" in c]]
        combo = self.frame[[c for c in self.frame if "freq" in c or "pool" in c]]
        if self.transform: covers = Clusterer.transforms[self.transform](covers)
        if self.transform: tetras = Clusterer.transforms[self.transform](tetras)
        if self.transform: combo = Clusterer.transforms[self.transform](combo)
        self.tetras_clusters = self.algorithm.fit_predict(tetras)        
        self.coverage_clusters = self.algorithm.fit_predict(covers)
        self.combo_clusters = self.algorithm.fit_predict(combo)
        self.clusters =  zip(self.frame.index,self.combo_clusters)
        if self.parent.clusters is None: self.parent.clusters = self.clusters
        #self.export()

    def export(self):
        parameters = "-".join(["nb_clusts",str(self.number_clusts),"min_len",str(self.min_length),"max_freq",str(self.max_freq),"transform",self.transform])
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
        path=os.path.join(self.p.combo,parameters)
        for i in set(self.combo_clusters):
            contig_list=[]
            for j in range(0,len(contigs)):
                if(self.combo_clusters[j] == i): contig_list.append(contigs[j])
            assembly.write_contiglist(contig_list,path,"bin_"+str(i)+".fasta")

