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
    Expects version 1.10"""

    short_name = "prokka"

    all_paths= """
    /output/
    /output/prokka_11112014.txt
    /output/prokka_11112014.ecn
    /output/prokka_11112014.err
    /output/prokka_11112014.faa
    /output/prokka_11112014.ffn
    /output/prokka_11112014.fna
    /output/prokka_11112014.fsa
    /output/prokka_11112014.gbf
    /output/prokka_11112014.gff
    /output/prokka_11112014.log
    /output/prokka_11112014.sqn
    /output/prokka_11112014.tbl
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
                  '--locustag', 'prokka',
                  '--compliant',
                  '--usegenus',
                  '--metagenome',
                  '--addgene',
                  '--quiet',
                  '--force',
                  self.contigs.fasta)

###############################################################################
class ProkkaResults(object):

    def __nonzero__(self): return self.prokka.p.log.exists
    def __init__(self, prokka):
        self.prokka = prokka

    @property_cached
    def proteins(self):
        """Returns ."""
        return dict(line.strip('\n').split(',') for line in self.concoct.p.clustering)
