# Built-in modules #
import os

# Internal modules #
from gefes.taxonomy.assignments import Assignment

# First party modules #
from plumbing.autopaths import AutoPaths, DirectoryPath, FilePath
from plumbing.cache     import property_cached
from plumbing.slurm     import num_processors
from plumbing.common    import which, SuppressAllOutput

# Third party modules #
import sh

# Prints warnings #
with SuppressAllOutput(): import ete2

###############################################################################
class Phylophlan(object):
    """Use Phylophlan to predict the taxonomy of bins. But beware:
    - Changelog stops at May 2013. Obsolete ?
    - Three versions available: Latest official release, latest tag or latest commit.
    - It requires usearch 5 to be in the PATH as `usearch` but doesn't check T_T
    - You have to manually change line 28 of the script after installation :'(
    - Cannot specify input and output directories... these are fixed, I kid you not.
    - Strangely changes behavior if no TTY is attached to its STDIN -.-
    - AttributeError in some cases when too few proteins inputted ?
    - You'd better create a directory named 'output' or it will crash (OSError)
    - All protein IDs have to be unique.
    - You have to remove the contents of the output directory before rerunning the tool.
    - The "low_conf" file has an underscore but not "medium-conf" and "high-conf".
    - Outputs windows line ends mixed with unix line ends on STDOUT.
    - Complains about '*' letter in FAA files.
    - Because of crazy input scheme only one instance can run at the same time."""

    short_name   = 'phylophlan'
    long_name    = 'PhyloPhlAn v1.1'
    executable   = 'phylophlan.py'
    url          = 'https://bitbucket.org/nsegata/phylophlan/'
    dependencies = ['muscle', 'usearch', 'FastTree']

    install_script = """
        cd ~/programs/
        hg clone https://bitbucket.org/nsegata/phylophlan
        cd phylophlan
        sed -i "s/usearch/usearch5/g" $HOME/programs/phylophlan/phylophlan.py
        escaped_home=${HOME//\//\\\/}
        sed -i "s/taxcuration\/')/$escaped_home\/programs\/phylophlan\/taxcuration\/')/g" $HOME/programs/phylophlan/phylophlan.py"""

    bash_rc = ["""export PATH="$HOME/programs/phylophlan":$PATH"""]

    all_paths = """
    /data/
    /input/
    /output/
    /stdout.txt
    /stderr.txt
    /output/pruned_tree.nwk
    """

    def __repr__(self): return '<%s object on binner %s>' % (self.__class__.__name__, self.binner)

    def __nonzero__(self):
        path = self.base_dir + 'output/' + self.proj_code + '/imputed_conf_low_conf.txt'
        return FilePath(path).exists

    def __init__(self, binner, result_dir):
        # Save attributes #
        self.binner, self.parent = binner, binner
        self.result_dir = result_dir
        # Auto paths #
        self.base_dir = self.result_dir + self.short_name + '/'
        self.p = AutoPaths(self.base_dir, self.all_paths)
        # Short project name #
        self.proj_code = self.binner.samples[0].project.name

    def run(self, cpus=None):
        # Variable threads #
        if cpus is None: cpus = num_processors
        # Where is the executable #
        program_dir = FilePath(which(self.executable)).directory
        # Link the data directory #
        self.p.data_dir.link_from(program_dir + 'data/', safe=True)
        # Crazy fixed input directory #
        input_dir = DirectoryPath(self.p.input_dir + self.proj_code)
        input_dir.remove()
        input_dir.create()
        # Link all input faa files #
        for b in self.binner.results.good_bins:
            b.faa.link_to(input_dir + "bin_" + b.name + '.faa')
        # Crazy fixed output directory #
        output_dir = self.p.output_dir
        output_dir.remove()
        output_dir.create()
        # Change directory #
        current_dir = os.getcwd()
        os.chdir(self.base_dir)
        # Call the executable #
        command = sh.Command("phylophlan.py")
        command('-i', # Integrates into the existing tree of life
                '-t', # Predicts taxonomy
                '--nproc', cpus,
                self.proj_code, # Name of the input directory
                _tty_in = True, # Without, it changes behavior
                _out = self.p.stderr.path,
                _err = self.p.stdout.path)
        # Make tree pruned tree #
        self.results.pruned_tree.save(outfile=self.p.pruned_tree)
        # Restore #
        os.chdir(current_dir)

    @property_cached
    def results(self):
        results = PhylophlanResults(self)
        if not results: raise Exception("You can't access results from Phylophlan before running the algorithm.")
        return results

###############################################################################
class PhylophlanResults(object):

    example_paths = """
    /output/proj/imputed_conf_low_conf.txt
    /output/proj/imputed_conf_high-conf.txt
    /output/proj/imputed_conf_medium-conf.txt
    /output/proj/incomplete_conf_high-conf.txt
    /output/proj/proj.tree.int.nwk"""

    def __nonzero__(self):
        path = self.base_dir + 'output/' + self.proj_code + '/imputed_conf_low_conf.txt'
        return FilePath(path).exists

    def __init__(self, phylophlan):
        self.phylophlan = phylophlan
        self.binner     = phylophlan.parent
        self.base_dir   = self.phylophlan.base_dir
        self.proj_code  = self.phylophlan.proj_code
        self.p          = self.phylophlan.p
        # The results directory #
        self.results_dir = DirectoryPath(self.base_dir + 'output/' + self.proj_code + '/')

    @property_cached
    def files(self):
        files = [self.base_dir + 'output/' + self.proj_code + '/imputed_conf_low_conf.txt',
                 self.base_dir + 'output/' + self.proj_code + '/imputed_conf_high-conf.txt',
                 self.base_dir + 'output/' + self.proj_code + '/imputed_conf_medium-conf.txt',
                 self.base_dir + 'output/' + self.proj_code + '/incomplete_conf_high-conf.txt',
                 self.base_dir + 'output/' + self.proj_code + '/incomplete_conf_medium-conf.txt',
                 self.base_dir + 'output/' + self.proj_code + '/incomplete_conf_low-conf.txt']
        files = map(FilePath, files)
        files = [f for f in files if f.exists]
        return files

    @property_cached
    def assignments(self):
        files = self.results_dir.glob('*.txt')
        lines = (line for f in files for line in f)
        def line_to_kv(line):
            bin_id, assignment = line.split('\t')[0:2]
            return bin_id, Assignment(assignment)
        return dict(map(line_to_kv, lines))

    @property
    def tree(self):
        return self.base_dir + 'output/' + self.proj_code + '/' + self.proj_code + '.tree.int.nwk'

    @property
    def tree_ete(self):
        """The tree as an object in python memory from ETE2."""
        return ete2.Tree(self.tree)

    @property
    def pruned_tree(self):
        """A tree containing every one of the good bins placed within
        the tree of life. We will build it by starting with a full tree
        of life and then pruning any node that is not connected to one
        of the bins."""
        # Prune it #
        tree   = self.tree_ete
        pruned = tree.prune(["bin_" + b.name for b in self.binner.results.good_bins], preserve_branch_length=True)
        pruned.save(outfile=self.p.pruned_tree)
        return pruned