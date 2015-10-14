# Built-in modules #
import os

# Internal modules #

# First party modules #
from plumbing.autopaths import AutoPaths, DirectoryPath, FilePath
from plumbing.cache     import property_cached
from plumbing.slurm     import num_processors
from plumbing.common    import which

# Third party modules #
import sh

###############################################################################
class Phylophlan(object):
    """Use Phylophlan to predict the taxonomy of bins.
    - Changelog stops at May 2013
    - It requires usearch version 5 to be in the PATH as `usearch` T_T
    - You have to manually change line 28 of the script after installation :'(
    - Has a different behavior if no TTY is attached to its STDIN :x"""

    short_name = 'phylophlan'
    long_name  = 'PhyloPhlAn v0.99'
    executable = 'phylophlan.py'
    url        = 'https://bitbucket.org/nsegata/phylophlan/'
    dependencies = ['muscle', 'usearch', 'FastTree']

    all_paths = """
    /data/
    /input/proj/proj.faa
    /output/proj/
    """

    def __repr__(self): return '<%s object on %s>' % (self.__class__.__name__, self.bin)

    def __init__(self, bin, result_dir):
        # Save attributes #
        self.bin = bin
        self.result_dir = result_dir
        # Auto paths #
        self.base_dir = self.result_dir + self.short_name + '/'
        self.p = AutoPaths(self.base_dir, self.all_paths)

    def run(self, cpus=None):
        # Variable threads #
        if cpus is None: cpus = num_processors
        # Crazy fixed input and output directories #
        working_dir = DirectoryPath(self.base_dir)
        working_dir.create(safe=True)
        program_dir = FilePath(which(self.executable)).directory
        self.p.proj_faa.link_from(self.bin.faa)
        self.p.data_dir.link_from(program_dir + 'data/', safe=True)
        current_dir = os.getcwd()
        os.chdir(working_dir)
        # Call the executable #
        command = sh.Command("phylophlan.py")
        command('--nproc', cpus, 'proj', _tty_in=True)
        # Restore #
        os.chdir(current_dir)

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