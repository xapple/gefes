# Built-in modules #

# Internal modules #

# First party modules #
from plumbing.autopaths import AutoPaths
from plumbing.slurm import num_processors
from plumbing.cache import property_cached

# Third party modules #
import sh

###############################################################################
class Prodigal(object):
    """Will run the Prodigal software on one single contig.
    Expects version  2.6.2
    Good read at https://github.com/hyattpd/prodigal/wiki/Advice-by-Input-Type"""

    short_name = 'prodigal'
    long_name  = 'Prodigal v2.6.2'
    executable = 'prodigal'
    url        = 'https://github.com/hyattpd/Prodigal'
    dependencies = []

    all_paths= """
    /output/coords.gbk
    /output/translated.faa
    """

    def __nonzero__(self): pass

    def __init__(self, contig, result_dir):
        # Save Attributes #
        self.contig = contig
        self.result_dir = result_dir
        # Auto paths #
        self.base_dir = self.result_dir + self.short_name + '/'
        self.p = AutoPaths(self.base_dir, self.all_paths)

    def run(self, cpus=None):
        # Variable threads #
        if cpus is None: cpus = num_processors
        # Run it #
        sh.prodigal('-i', self.contig.fasta,
                   '-o', self.p.gbk,
                   '-a', self.p.faa,
                   '-p', 'anon')

    @property_cached
    def results(self):
        results = ProdigalResults(self)
        if not results: raise Exception("You can't access results from Prodigal before running the algorithm.")
        return results

###############################################################################
class ProdigalResults(object):

    def __nonzero__(self): pass
    def __init__(self, prodigal):
        self.prodigal = prodigal
