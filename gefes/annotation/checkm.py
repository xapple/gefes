# Built-in modules #

# Internal modules #
import sys

# First party modules #
from plumbing.autopaths import AutoPaths
from plumbing.cache import property_cached
from plumbing.slurm import num_processors

# Third party modules #
import sh

###############################################################################
class Checkm(object):
    """Use CheckM at https://github.com/Ecogenomics/CheckM/wiki
    to evaluate a bin of contigs.
    Expects version v0.9.7.
    """

    short_name = 'checkm'
    long_name  = 'CheckM v0.9.7'
    executable = 'checkm'
    dependencies = ['hmmer', 'prodigal', 'pplacer']

    all_paths = """
    /contigs.fasta
    /stdout.txt
    /stderr.txt
    /output/
    """

    def __init__(self, bin, result_dir):
        # Save attributes #
        self.bin = bin
        # Auto paths #
        self.result_dir = result_dir
        self.p = AutoPaths(self.result_dir, self.all_paths)

    def run(self):
        # Link the bin's fasta file #
        self.bin.fasta.link_to(self.p.fasta)
        # Run the pipeline #
        print "Launching CheckM on bin '%s'..." % self.bin.name; sys.stdout.flush()
        sh.checkm('lineage_wf',
                  '-x', 'fasta',
                  '-t', num_processors,
                  self.result_dir,
                  self.p.output_dir,
                  _out=self.p.stdout.path,
                  _err=self.p.stderr.path)

    @property_cached
    def results(self):
        results = CheckmResults(self)
        if not results: raise Exception("You can't access results from CheckM before running the algorithm.")
        return results

###############################################################################
class CheckmResults(object):

    def __nonzero__(self): return bool(self.checkm.p.stdout)

    def __init__(self, checkm):
        self.checkm = checkm

    @property_cached
    def statistcs(self):
        """The various statistics produced by checkm in a dictionary."""
        keys = ["bin_id", "lineage", "genomes", "markers", "marker_sets",
                "0", "1", "2", "3", "4", "5+",
                "completeness", "contamination", "heterogeneity"]
        return dict(zip(keys, list(self.checkm.p.stdout)[4].split()))