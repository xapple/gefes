# Built-in modules #
import warnings

# Internal modules #

# First party modules #
from plumbing.autopaths import AutoPaths
from plumbing.cache import property_cached
from plumbing.slurm import num_processors
from seqsearch.databases.pfam import pfam

# Warnings #
warnings.simplefilter("ignore", "Bio.SearchIO")
warnings.simplefilter("ignore", "BiopythonWarning")
from Bio import SearchIO

# Third party modules #
import sh

###############################################################################
class Hmmer(object):
    """Takes care of running Hmmer."""

    short_name = 'hmmer'
    long_name  = 'HMMER 3.1b2 (February 2015)'
    executable = 'hmmsearch'
    url        = 'http://hmmer.org/'
    license    = 'GPLv3'
    dependencies = []

    all_paths= """
    /seq_hits.txt
    """

    def __nonzero__(self): return bool(self.p.hits)

    def __init__(self, proteins, result_dir, database='pfam', e_value=10**-5):
        # Save Attributes #
        self.proteins   = proteins
        self.result_dir = result_dir
        self.database   = database
        # Thresholds #
        self.e_value = 10**-5
        # Auto detect database #
        if self.database == 'pfam': self.database = pfam.hmm_db
        # Auto paths #
        self.base_dir = self.result_dir + self.short_name + '/'
        self.p = AutoPaths(self.base_dir, self.all_paths)

    @property
    def command_args(self):
            return (
                '-o', '/dev/null',
                '--tblout', self.p.hits, # parseable table of per-sequence hits
                '--notextw',        # unlimited ASCII text output line width
                '--acc',            # prefer accessions over names in output
                '--seed', 1,        # set RNG seed to <n>
                '-E', self.e_value, # report only sequences <= this e-value
                self.database,
                self.proteins,
            )

    def run(self, cpus=None):
        # Check if FASTA is not empty #
        if self.proteins.count_bytes == 0:
            warnings.warn("Hmmer search on a file with no proteins", RuntimeWarning)
            return False
        # Variable threads #
        if cpus is None: cpus = num_processors
        # Run it #
        sh.hmmsearch('--cpu', cpus, *self.command_args)

    @property_cached
    def results(self):
        results = HmmerResults(self)
        if not results: raise Exception("You can't access results from Hmmer before running the algorithm.")
        return results

###############################################################################
class HmmerResults(object):

    def __nonzero__(self): return bool(self.hmmer.p.hits)

    def __init__(self, hmmer):
        self.hmmer = hmmer
        self.p     = hmmer.p

    @property
    def hits(self):
        return SearchIO.read(self.p.hits, 'hmmer3-tab')
