# Futures #
from __future__ import division

# Built-in modules #

# Internal modules #
from gefes.common.autopaths import AutoPaths
from gefes.common.cache import property_cached

# Third party modules #
import pandas
from sklearn import cluster
from scipy.spatial import distance
from scipy.cluster import hierarchy
from sklearn.cluster  import KMeans

###############################################################################
class Clusterer(object):
    """Recieves the matrix detailing all information for every contig,
    and is responsible for deciding which contigs go together."""

    all_paths = """
    /lorem.txt
    """

    def __repr__(self): return '<%s object of %s>' % (self.__class__.__name__, self.parent)

    def __init__(self, pool):
        # Save parent #
        self.parent, self.pool = pool, pool
        # Auto paths #
        self.base_dir = self.parent.p.clustering
        self.p = AutoPaths(self.base_dir, self.all_paths)
        self.kmeans = GefesKMeans(self)
        self.min_length = 0;
        self.max_freq = 0;
        
    def run(self):
        pass

    def frame_filter(self):
        temp_frame = self.parent.frame
        if(self.min_length):
            temp_frame = temp_frame[temp_frame.length > self.min_length]
        if(self.max_freq):
            good_ones = temp_frame[[c for c in temp_frame if "freq" in c]].apply(lambda x: sum(x>self.max_freq)!=0,1)
            temp_frame = temp_frame[good_ones]
        tetras = temp_frame[[c for c in temp_frame if "freq" in c]]
        covers = temp_frame[[c for c in temp_frame if "freq" not in c and c!="length"]]
        names = [f[0] for f in temp_frame.itertuples()]
        return (names,tetras,covers)

    
###############################################################################
class GefesKMeans(object):
    """Receives the matrix and uses kmeans to cluster it in N different (linearly separated) clusters"""

    def __init__(self,parent,number_clusts=8):
        self.number_clusts=number_clusts
        self.parent = parent

    @property_cached
    def algorithm(self):
        return KMeans(self.number_clusts)

    @property_cached
    def tetras_clusters(self):
        (names,tetras,covers) = self.parent.frame_filter()
        return self.algorithm.fit_predict(tetras)

    @property_cached
    def coverage_clusters(self):
        (names,tetras,covers) = self.parent.frame_filter()
        return self.algorithm.fit_predict(covers)

    
