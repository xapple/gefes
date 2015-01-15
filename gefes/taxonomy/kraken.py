# Built-in modules #
import os

# Internal modules #

# First party modules #
from plumbing.autopaths import AutoPaths
from plumbing.cache import property_cached
from plumbing.slurm import num_processors

# Third party modules #
import sh

# Constants #
home = os.environ['HOME'] + '/'

###############################################################################
class Kraken(object):
    """Use Kraken at http://bowtie-bio.sourceforge.net/bowtie2/index.shtml
    to predict taxonomy on the raw reads.
    """

    all_paths = """
    /output.txt
    """

    def __repr__(self): return '<%s object on %s>' % (self.__class__.__name__, self.parent)

    def __init__(self, reads):
        # Save attributes #
        self.reads = reads
        # Auto paths #
        self.base_dir = self.result_dir + 'bowtie/'
        self.p = AutoPaths(self.base_dir, self.all_paths)

    def run(self):
        sh.kraken('--preload', '--fastq-input', '--gzip-compressed',
                  '--threads', num_processors,
                  '--db',      home + 'databases/kraken/standard',
                  '--output',  self.p.output,
                  self.reads)

    @property_cached
    def results(self):
        results = KrakenResults(self)
        if not results: raise Exception("You can't access results from Kraken before running the tool.")
        return results

###############################################################################
class KrakenResults(object):

    all_paths = """
    /output/lorem
    """

    def __nonzero__(self): return 0
    def __init__(self, kraken):
        self.kraken = kraken
        self.p = AutoPaths(self.kraken.base_dir, self.all_paths)

    @property_cached
    def composition(self):
        return self.p.output.content