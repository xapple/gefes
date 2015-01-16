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
standard_db = home + 'databases/kraken/standard'

###############################################################################
class Kraken(object):
    """Use Kraken at http://bowtie-bio.sourceforge.net/bowtie2/index.shtml
    to predict taxonomy on the raw reads.
    """

    all_paths = """
    /raw_output.txt
    /summary.tsv
    """

    def __repr__(self): return '<%s object on %s>' % (self.__class__.__name__, self.parent)

    def __init__(self, source, base_dir=None):
        # Basic #
        self.source   = source
        self.base_dir = base_dir
        # Default case #
        if base_dir is None: self.base_dir = self.source.fwd.prefix_path + '.kraken/'
        # Auto paths #
        self.p = AutoPaths(self.base_dir, self.all_paths)

    def run(self, keep_raw=False):
        # Run the main classification #
        sh.kraken('--preload', '--fastq-input', '--gzip-compressed',
                  '--threads', str(num_processors),
                  '--db',      standard_db,
                  '--output',  self.p.raw_output,
                  '--paired',  self.source.fwd, self.source.rev)
        # Run the report #
        report = sh.Command('kraken-report')
        report('--db', standard_db, self.p.raw_output, _out=str(self.p.summary))
        # Remove intermediary output #
        if not keep_raw: self.p.raw_output.remove()

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

    def __nonzero__(self): return bool(self.p.summary)
    def __init__(self, kraken):
        self.kraken = kraken
        self.p = AutoPaths(self.kraken.base_dir, self.all_paths)

    @property_cached
    def composition(self):
        return self.p.output.content

