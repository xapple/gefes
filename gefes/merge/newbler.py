# Built-in modules #
import sys

# Internal modules #
from gefes.merge import Merger
from gefes.assemble.contig import Contig

# First party modules #
from plumbing.autopaths import AutoPaths
from plumbing.cache import property_cached

# Third party modules #
import sh
from shell_command import shell_output

# Constant #

###############################################################################
class Newbler(Merger):
    """Used for merging of several assemblies together.
    Acts as a new assembly.
    Newbler is outdated as Roche has stopped making sequencing machines in 2013.
    Newbler is closed source and Roche dismissed the community
    petition to make it open source.
    You can only obtain the software via a registration form."""

    short_name = 'newbler'
    long_name  = 'Newbler in Roche GS De Novo Assembler software v2.9'
    executable = 'runAssembly'
    dependencies = []

    all_paths = """
    /combined_cut_up.fasta
    /output/
    """

    def __init__(self, assemblies):
        # Base parameters #
        pass
        # Auto paths #
        self.base_dir = self.result_dir + self.short_name + '/' + str(self.kmer_size) + '/'
        self.p = AutoPaths(self.base_dir, self.all_paths)

    def run(self, verbose=True):
        # Check that all the sets of contigs have a cut-up version #
        for assembly in self.assemblies:
            if not assembly.p.cut_up:
                if verbose: print 'Cutting up contigs from %s' % assembly; sys.stdout.flush()
                sh.cut_up_fasta(assembly.results.contigs_fasta.path, _out=assembly.p.cut_up.path)
        # Merge the cut-uped fastas #
        if verbose: print 'Combining contigs from %s' % self.assemblies; sys.stdout.flush()
        paths = [assembly.p.cut_up.path for assembly in self.assemblies]
        shell_output('cat %s > %s' % (' '.join(paths), self.p.combined))
        # Call newbler #
        if verbose: print 'Calling Newbler...'; sys.stdout.flush()
        sh.runAssembly('-force', '-o', self.p.output_dir.path, self.p.combined.path)

    @property_cached
    def results(self):
        results = NewblerResults(self)
        if not results: raise Exception("You can't access results from Newbler before running the algorithm.")
        return results

###############################################################################
class NewblerResults(object):

    def __nonzero__(self): return 0 # bool(self.contigs_fasta)
    def __init__(self, newbler):
        self.newbler = newbler

    @property_cached
    def contigs(self):
        """All the contigs produced returned as a list of our Contig custom objects."""
        return [Contig(self.ray, record, num=i) for i,record in enumerate(self.contigs_fasta)]

    @property_cached
    def contig_id_to_contig(self):
        """A dictionary with contig names as keys and contig objects as values."""
        return {c.name: c for c in self.contigs}