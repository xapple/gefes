# Built-in modules #

# Internal modules #

# First party modules #
from plumbing.autopaths import AutoPaths
from plumbing.slurm import nr_threads

# Third party modules #
import sh

###############################################################################
class Prokka(object):
    """Will run the Prokka software on a set of contigs.
    http://www.vicbioinformatics.com/software.prokka.shtml
    Expects version 1.7.2"""

    short_name = "prokka"

    all_paths= """
    /output/
    """

    def __init__(self, contigs, result_dir):
        # Save Attributes #
        self.contigs = contigs
        self.result_dir = result_dir
        # Auto paths #
        self.base_dir = self.result_dir + self.short_name + '/'
        self.p = AutoPaths(self.base_dir, self.all_paths)

    def run(self):
        sh.prokka('--outdir', self.p.output,
                  '--cpus', nr_threads,
                  '--compliant',
                  '--usegenus',
                  '--metagenome',
                  '--addgene',
                  '--quiet',
                  '--force',
                  self.contigs.fasta)
