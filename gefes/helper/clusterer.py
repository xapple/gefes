# Futures #
from __future__ import division

# Built-in modules #

# Internal modules #
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
    """Recieves an assembly, outputs a clustering"""

    def __repr__(self): return '<%s object of %s>' % (self.__class__.__name__, self.parent)

    def __init__(self):
        self.clusters = None
                
    def log10(x):
        return x.applymap(math.log10)

    def log10p1(x):
        log10p = lambda z: math.log10(z+1)
        return x.applymap(log10p)

    def rank(x):
        return x.apply(scipy.stats.rankdata,0)

    transforms = {'log10' : log10, 'log10p1' : log10p1 , 'rank' : rank}
    
        
###############################################################################
class GefesKMeans(Clusterer):
    """Receives the matrix and uses kmeans to cluster it in N different (linearly separated) clusters"""
    
    def __init__(self,nb = 8,method = 'tetramer', max_freq = None,min_length = None, transform = None):
        # Save parent
        super(Clusterer,self).__init__()
        # Kmeans and params
        self.number_clusts=nb
        self.algorithm = KMeans(self.number_clusts)
        # Filters
        self.method = method
        self.max_freq = max_freq
        self.min_length = min_length
        self.transform = transform

    def run(self,assembly):
        self.frame = assembly.filtered_frame(self.max_freq,self.min_length)
        if self.method == 'tetramer':
            data = self.frame[[c for c in self.frame if "freq" in c]]
        elif self.method == 'coverage':
            data = self.frame[[c for c in self.frame if "pool" in c]]
        elif self.method == 'uwcombo':
            data = self.frame[[c for c in self.frame if "freq" in c or "pool" in c]]
        else: raise Exception("This method does not exist in da Kmeans")
                
        if self.transform: data = Clusterer.transforms[self.transform](data)

        self.clusters = zip(self.frame.index,self.algorithm.fit_predict(data))

