# Futures #
from __future__ import division

# Built-in modules #
from math import log10
import os
import json

# Internal modules #
from gefes.common.autopaths import AutoPaths
from gefes.helper.clusterer import Clusterer
import gefes.helper.clusterer
from gefes.common.cache import property_cached
from gefes.fasta.single import FASTA
from gefes.helper.contig import Contig
from gefes.helper.bin import Bin
from gefes.common.autopaths import AutoPaths
from gefes.running import Runner
from gefes.common.slurm import SLURMJob
from gefes.helper.clusterer import *
 
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
    def __iter__(self): return iter(self.binnings.values())
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
        self.binnings={}
        binnings_paths = os.listdir(self.base_dir)
        for b in binnings_paths:
            self.binnings[b]=Binning(self, b)
        


    def load(self):
        for b in self: b.load()
            
    def new(self,name,clusterer = {'type' : 'GefesKMeans', 'args' : {'nb' : 8 ,
                                                                    'method' : 'tetramer',
                                                                    'max_freq' : 0.1 ,
                                                                    'min_length' : 1000 ,
                                                                    'transform' : 'rank'
                                                                    }}):
        self.binnings[name] = Binning(self, name, clusterer = clusterer)
             

###############################################################################        
class Binning(object):

    all_paths = """
    /bins/
    /settings.json
    /linkage_matrix.csv
    /coverage_stats.csv
    /frame.csv
    /logs/
    """
    def __repr__(self): return '<'+self.__class__.__name__+ ' object with ' + str(len(self)) + ' bins>' 
    def __iter__(self): return iter(self.bins.values())
    def __len__(self):
        if hasattr(self, 'bins'):
             return len(self.bins)
        else:
            return 0
    def __getitem__(self, key):
        if isinstance(key,int):
            return self.bins.values()[key]
        else:
            return self.bins[key]
    
    def __init__(self, parent,name,clusterer=None):
        # Save parent #
        self.parent =  parent
        self.name = name
        # Auto paths #
        self.base_dir = self.parent.p._base_dir+name
        self.p = AutoPaths(self.base_dir, self.all_paths)
        # Runner #
        self.runner=BinningRunner(self)
        # Clusterer representation #
        if clusterer is not None:
            with open(self.p.settings, 'w') as outfile:
                    json.dump(clusterer, outfile)
        self.loaded = False 

    def load(self):
        self.frame = DataFrame.from_csv(self.p.frame)
        files = os.listdir(self.p.bins)
        self.bins={}
        for f in files:
            bini=Bin.fromfolder(self,f)
            self.bins[bini.name]=bini
        self.loaded =True

                            
    def cluster(self):
        self.clusterer = getattr(gefes.helper.clusterer, self.clusterer_rep['type'])(self.clusterer_rep['args'])
        self.clusterer.run(self.parent.assembly)
        self.bins = {}
        for contig_name,cluster in self.clusterer.clusters:
            cluster=str(cluster)
            contig=[o for o in self.parent.assembly.contigs if o.name is contig_name]
            if self.bins.has_key(cluster):
                    self.bins[cluster].extend(contig)
            else:
                self.bins[cluster]=Bin(self,contig,cluster)
        self.frame = self.clusterer.frame
        self.export()
        self.loaded = True

    def annotate(self):
        if not self.loaded: self.load()
        for bini in self:
            bini.annotate()
            
    @property_cached            
    def clusterer_rep(self):
        with open(self.p.settings) as j:
            clusterer_rep = json.load(j)
        return clusterer_rep

            
    @property_cached
    def linkage_matrix(self):
        if os.path.exists(self.p.linkage_matrix):
            return DataFrame.from_csv(self.p.linkage_matrix)
        else :
            linkage=[p.mapper.linkage for p in self.parent.parent]
        
            matrix = DataFrame(0,index=[b.name for b in self],columns=[b.name for b in self])
            for b1 in self:
                for b2 in self:
                    for contig1 in b1:
                        for contig2 in b2:
                            temp = 0
                            for p in linkage:
                                if contig1.name != contig2.name:
                                    matrix[b1.name][b2.name] = matrix[b1.name][b2.name] + sum(p[contig1.name][contig2.name])/2.0
                                    matrix[b2.name][b1.name] = matrix[b2.name][b1.name] + sum(p[contig1.name][contig2.name])/2.0
            return matrix

    @property_cached
    def coverage_stats(self):
        if os.path.exists(self.p.coverage):
            return DataFrame.from_csv(self.p.coverage)
        else :
            out={}
            for b in self:
                out[b.name]=self.frame[[col for col in self.frame if "pool" in col]].loc[[c.name for c in b]]
                out[b.name]=out[b.name].median()
        return DataFrame.from_dict(out)


    def export(self):
        for b in self:
            b.export()
        self.coverage_stats.to_csv(self.p.coverage)
        self.linkage_matrix.to_csv(self.p.linkage_matrix)
        self.frame.to_csv(self.p.frame)

###############################################################################

class BinningRunner(Runner):
    """Will run stuff on a project"""
    default_time = '7-00:00:00'

    default_steps = [
        {'cluster':    {}},
        {'annotate':   {}},
    ]

    def __init__(self, parent):
        # Save parent #
        self.parent, self.binning = parent, parent
        self.project = self.binning.parent.parent

    def run_slurm(self, steps=None, **kwargs):
        # Make script #
        command = """steps = %s
                     binning = [binning for binning in gefes.projects['%s'].binner if binning.name=='%s'][0]

                     binning.runner(steps)""" % (steps,self.project.name, self.binning.name)
        # Test case #
        if 'test' in self.project.name:
            kwargs['time'] = '00:15:00'
            kwargs['qos'] = False
            kwargs['email'] = '/dev/null'

        # Send it #
        if 'time' not in kwargs: kwargs['time'] = self.default_time
        if 'email' not in kwargs: kwargs['email'] = None
        job_name = "gefes_%s_%s" % (self.project.name,self.binning.name)
        self.slurm_job = SLURMJob(command, self.binning.p.logs_dir, job_name=job_name, **kwargs)
        self.slurm_job.launch()
