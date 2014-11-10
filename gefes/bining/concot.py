# Built-in modules #

# Internal modules #

# First party modules #
from plumbing.autopaths import AutoPaths
from plumbing.cache import property_cached

# Third party modules #

###############################################################################
class Concoct(object):
    """Use CONCOCT at https://github.com/BinPro/CONCOCT
    to bin contigs togather.
    """

    all_paths = """
    /lorem.fasta
    """

    def __repr__(self): return '<%s object on %s>' % (self.__class__.__name__, self.parent)

    def __init__(self, contigs):
        # Save attributes #
        self.contigs = contigs
        # Auto paths #
        self.base_dir = self.result_dir + 'bowtie/'
        self.p = AutoPaths(self.base_dir, self.all_paths)

    def run(self):
        pass

    @property_cached
    def results(self):
        results = ConcoctResults(self)
        if not results: raise Exception("You can't access results from ConcoctResults before running the binning.")
        return results

###############################################################################
class ConcoctResults(object):

    all_paths = """
    /output/lorem
    """

    def __nonzero__(self): return 0
    def __init__(self, concoct):
        self.concoct = concoct
        self.p = AutoPaths(self.concoct.base_dir, self.all_paths)

    @property_cached
    def lorem(self):
        pass