# Built-in modules #

# Internal modules #

# First party modules #
from plumbing.autopaths import AutoPaths
from plumbing.cache import property_cached
from fasta import FASTA

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
    /output/proteins.fna
    /output/translated.faa
    """

    def __nonzero__(self): return bool(self.p.proteins)

    def __init__(self, contig, result_dir):
        # Save Attributes #
        self.contig = contig
        self.result_dir = result_dir
        # Auto paths #
        self.base_dir = self.result_dir + self.short_name + '/'
        self.p = AutoPaths(self.base_dir, self.all_paths)

    @property
    def command_args(self):
        return ('-i', self.contig.fasta,
                '-o', self.p.gbk,
                '-d', self.p.fna,
                '-a', self.p.faa,
                '-p', 'meta')

    def run(self):
        sh.prodigal(*self.command_args)

    @property_cached
    def results(self):
        results = ProdigalResults(self)
        if not results: raise Exception("You can't access results from Prodigal before running the algorithm.")
        return results

###############################################################################
class ProdigalResults(object):

    def __nonzero__(self): return bool(self.prodigal.p.proteins)
    def __iter__(self): return iter(self.faa)

    def __init__(self, prodigal):
        self.prodigal = prodigal

    @property_cached
    def fasta(self): return FASTA(self.prodigal.p.fna)

    @property_cached
    def faa(self): return FASTA(self.prodigal.p.faa)
