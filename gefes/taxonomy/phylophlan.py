# Built-in modules #

# Internal modules #

# First party modules #
from plumbing.autopaths import AutoPaths
from plumbing.cache import property_cached
from plumbing.slurm import num_processors

# Third party modules #
import sh

###############################################################################
class Phylophlan(object):
    """Use Phylophlan to predict the taxonomy of bins.
    - Changelog stops at May 2013
    - It requires usearch v5 to be in the path as `usearch` T_T
    - You have to manually change line 28 of the script T_T"""

    short_name = 'phylophlan'
    long_name  = 'PhyloPhlAn v0.99'
    executable = 'phylophlan.py'
    url        = 'https://bitbucket.org/nsegata/phylophlan/'
    dependencies = ['muscle', 'usearch', 'FastTree']

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

    def run(self, cpus=None):
        # Variable threads #
        if cpus is None: cpus = num_processors
        # Call the executable #
        command = sh.Command(self.executable)
        command('a', '--nproc', cpus)

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