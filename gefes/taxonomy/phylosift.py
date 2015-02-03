# Built-in modules #

# Internal modules #

# First party modules #
from plumbing.autopaths import AutoPaths
from plumbing.cache import property_cached

# Third party modules #
import sh

###############################################################################
class Phylosift(object):
    """Use Phylosift at https://github.com/gjospin/PhyloSift
    to conduct phylogenetic analysis.
    """

    short_name = 'phylosift'
    long_name  = 'Phylosift v20140419'
    executable = 'phylosift'
    dependencies = []

    all_paths = """
    /results/
    """

    def __repr__(self): return '<%s object on %s>' % (self.__class__.__name__, self.contig)

    def __init__(self, contig, result_dir):
        # Save attributes #
        self.contig = contig
        self.result_dir = result_dir
        # Auto paths #
        self.base_dir = self.result_dir + self.short_name + '/'
        self.p = AutoPaths(self.base_dir, self.all_paths)

    def run(self):
        sh.phylosift('all', '--output=' + self.base_dir, self.contig.fasta)

    @property_cached
    def results(self):
        results = PhylosiftResults(self)
        if not results: raise Exception("You can't access results from Phylosift before running the tool.")
        return results

###############################################################################
class PhylosiftResults(object):

    all_paths = """
    /output/lorem
    """

    def __nonzero__(self): return 0
    def __init__(self, phylosift):
        self.phylosift = phylosift
        self.p = AutoPaths(self.phylosift.base_dir, self.all_paths)

    @property_cached
    def lorem(self):
        pass