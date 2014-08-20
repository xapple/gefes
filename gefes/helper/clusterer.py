# Futures #
from __future__ import division

# Built-in modules #

# Internal modules #
from fasta import FASTA
from plumbing.autopaths import AutoPaths

# Third party modules #
import os
from pandas import DataFrame
import scipy.stats
from sklearn.cluster  import KMeans
import math

###############################################################################
class Clusterer(object):
    """Receives an assembly, outputs a clustering"""

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

    def __init__(self,parent, args = {'nb' : 8, 'method' : 'tetramer', 'max_freq' : None,'min_length' : None, 'transform' : None} ):
        # Save parent
        super(Clusterer,self).__init__()
        # Kmeans and params
        self.number_clusts=args['nb']
        self.algorithm = KMeans(self.number_clusts)
        # Filters
        self.method = args['method']
        self.max_freq = args['max_freq']
        self.min_length = args['min_length']
        self.transform = args['transform']

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

###############################################################################
class GefesCONCOCT(Clusterer):

    all_paths = """
    /clustering_gt1000.csv
    /coverages.csv
    /contigs.fasta
    """

    def __init__(self,parent, args = {'max_clusters' : 100, 'max_freq' : None, 'min_length' : None, 'transform' : None} ):
        # Save parent
        super(Clusterer,self).__init__()
        self.parent = parent
        # Auto paths #
        self.base_dir = self.parent.p._base_dir + "/concoct/"
        self.p = AutoPaths(self.base_dir, self.all_paths)
        # Filters
        self.min_length = args['min_length']
        if args.has_key('transform') : self.transform = args['transform']
        else : self.transform = None
        if args.has_key('max_clusters') : self.max_clusters = args['max_clusters']
        else :  self.max_clusters = 100
        if args.has_key('max_freq') : self.max_freq = args['max_freq']
        else : self.max_freq=1

    def run(self,assembly):
        self.frame = assembly.filtered_frame(self.max_freq,self.min_length)
        data = self.frame[[c for c in self.frame if "pool" in c]]
        if self.transform: data = Clusterer.transforms[self.transform](data)
        data.to_csv(self.p.coverages,sep="\t")
        self.contigs = data.index
        with FASTA(self.p.contigs) as ffile:
            ffile.write([seq.record for seq in assembly.contigs if seq.name in self.contigs])
        cwd = os.getcwd()
        os.chdir(self.p._base_dir)
        os.system("concoct " + " --coverage_file " + self.p.coverages + " --composition_file " +  self.p.contigs + " -c " + str(self.max_clusters))
        os.chdir(cwd)
        tmp_data = DataFrame.from_csv(self.p.clustering,header=None)
        self.clusters = zip(tmp_data.index,[int(v) for v in tmp_data[[1]].values])


###############################################################################

class GefesImport(Clusterer):
    """Importing an externally clustering with a simple tab separated file with the first column the contig names, the second the clusters (numbers from 0 to N) """

    all_paths = """
    /inport.csv
    """

    def __init__(self,parent, args = {'file' : "import.csv"} ):
        # Save parent
        super(Clusterer,self).__init__()
        self.parent = parent
        # Auto paths #
        self.base_dir = self.parent.p._base_dir + "/import/"
        self.p = AutoPaths(self.base_dir, self.all_paths)
        # Filters
        self.file = args['file']

    def run(self,assembly):
        self.frame = assembly.filtered_frame(1,0)
        os.system("cp " + self.file + " " + self.p.inport)
        tmp_data = DataFrame.from_csv(self.p.inport,header=None,sep="\t")
        self.contigs = tmp_data.index
        self.frame=self.frame.loc[list(tmp_data.index)]
        self.clusters = zip(tmp_data.index,[int(v) for v in tmp_data[[1]].values])

