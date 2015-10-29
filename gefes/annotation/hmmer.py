# Built-in modules #

# Internal modules #

# First party modules #
from plumbing.autopaths import AutoPaths
from plumbing.cache import property_cached
from plumbing.slurm import num_processors
from seqsearch.databases.pfam import pfam

# Third party modules #
import sh

###############################################################################
class Hmmer(object):
    """Takes care of running Hmmer"""

    short_name = 'hmmer'
    long_name  = 'HMMER 3.1b2 (February 2015)'
    executable = 'hmmsearch'
    url        = 'http://hmmer.org/'
    license    = 'GPLv3'
    dependencies = []

    all_paths= """
    /lorem
    """

    def __nonzero__(self): return self.p.proteins.exists

    def __init__(self, proteins, result_dir, database='pfam'):
        # Save Attributes #
        self.proteins   = proteins
        self.result_dir = result_dir
        self.database   = database
        # Auto detect database #
        if self.database == 'pfam': self.database = pfam.hmm_db
        # Auto paths #
        self.base_dir = self.result_dir + self.short_name + '/'
        self.p = AutoPaths(self.base_dir, self.all_paths)

    @property
    def command_args(self):
        return ('-o', self.p.gbk,
                '--notextw', # unlimit ASCII text output line width
                '--acc',     # prefer accessions over names in output
                '--seed', 1, # set RNG seed to <n>
                self.database,
                self.proteins)

    def run(self, cpus=None):
        # Variable threads #
        if cpus is None: cpus = num_processors
        # Run it #
        sh.hmmer('--cpu', cpus, *self.command_args)

    @property_cached
    def results(self):
        results = HmmerResults(self)
        if not results: raise Exception("You can't access results from Hmmer before running the algorithm.")
        return results

###############################################################################
class HmmerResults(object):

    def __nonzero__(self): return self.hmmer.p.proteins.exists
    def __iter__(self): return iter(self.faa)

    def __init__(self, hmmer):
        self.hmmer = hmmer
