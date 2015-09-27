# Built-in modules #

# Internal modules #

# First party modules #
from fasta import FASTA
from plumbing.autopaths import AutoPaths
from plumbing.slurm import num_processors
from plumbing.cache import property_cached

# Third party modules #
import sh

###############################################################################
class Prokka(object):
    """Will run the Prokka software on one single contig.
    Expects version 1.11"""

    short_name = 'prokka'
    long_name  = 'Prokka v1.11'
    executable = 'prokka'
    url        = 'http://www.vicbioinformatics.com/software.prokka.shtml'
    dependencies = ['parallel']

    all_paths= """
    /output/
    /output/L_*.txt
    /output/L_*.ecn
    /output/L_*.err
    /output/L_*.faa
    /output/L_*.ffn
    /output/L_*.fna
    /output/L_*.fsa
    /output/L_*.gbf
    /output/L_*.gff
    /output/L_*.log
    /output/L_*.sqn
    /output/L_*.tbl
    """

    def __nonzero__(self): return not self.p.output_dir.empty

    def __init__(self, contig, result_dir):
        # Save Attributes #
        self.contig = contig
        self.result_dir = result_dir
        # Auto paths #
        self.base_dir = self.result_dir + self.short_name + '/'
        self.p = AutoPaths(self.base_dir, self.all_paths)

    def run(self, cpus=None):
        """There is a subtle bug that can appear here when the fasta sequence
        header names are longer than a few characters.
        See https://github.com/tseemann/prokka/issues/135"""
        # Variable threads #
        if cpus is None: cpus = num_processors
        # Run it #
        sh.prokka('--outdir', self.p.output,
                  '--cpus', cpus,
                  '--locustag', 'L',
                  '--centre', 'C',
                  '--compliant',
                  '--usegenus',
                  '--metagenome',
                  '--addgene',
                  '--quiet',
                  '--force',
                  self.contig.fasta)

    @property_cached
    def results(self):
        results = ProkkaResults(self)
        if not results: raise Exception("You can't access results from Prokka before running the annotation.")
        return results

###############################################################################
class ProkkaResults(object):

    def __nonzero__(self): return not self.prokka.p.output_dir.empty
    def __init__(self, prokka):
        self.prokka = prokka

    @property_cached
    def functions(self):
        """Returns a list of predicted functions, one per predicted protein."""
        return [p.description[len(p.id)+1:] for p in FASTA(self.prokka.p.faa)]

    @property_cached
    def species(self):
        """Returns the predicted species of this contig."""
        return [p.description[len(p.id)+1:] for p in FASTA(self.prokka.p.fsa)][0]