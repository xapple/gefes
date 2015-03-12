# Built-in modules #

# Internal modules #

# First party modules #
from plumbing.autopaths import AutoPaths
from plumbing.cache import property_cached

# Third party modules #

###############################################################################
class Phylophlan(object):
    """Use Phylophlan at http://example.com
    to predict taxonomy on bins.
    """

    all_paths = """
    /lorem.fasta
    """

    def __repr__(self): return '<%s object on %s>' % (self.__class__.__name__, self.bin)

    def __init__(self, bin, result_dir):
        # Save attributes #
        self.bin = bin
        self.result_dir = result_dir
        # Auto paths #
        self.base_dir = self.result_dir + 'phylophlan/'
        self.p = AutoPaths(self.base_dir, self.all_paths)

    def run(self):
        pass

    @property_cached
    def results(self):
        results = PhylophlanResults(self)
        if not results: raise Exception("You can't access results from Phylophlan before running the algorithm.")
        return results

###############################################################################
class PhylophlanResults(object):

    all_paths = """
    /output/lorem
    """

    def __nonzero__(self): return 0
    def __init__(self, phylophlan):
        self.phylophlan = phylophlan
        self.p = AutoPaths(self.phylophlan.base_dir, self.all_paths)

    @property_cached
    def lorem(self):
        pass