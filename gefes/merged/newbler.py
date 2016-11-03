# Built-in modules #
import sys, os
from collections import OrderedDict

# Internal modules #
from gefes.merged import Merger
from gefes.assemble        import AssemblyResults
from gefes.report.assembly import AssemblyReport

# First party modules #
from plumbing.autopaths import AutoPaths
from plumbing.cache import property_cached
from plumbing.slurm import num_processors
from plumbing.common import andify
from fasta import FASTA

# Third party modules #
import sh
from shell_command import shell_output

# Constant #

###############################################################################
class Newbler(Merger):
    """Used for merging of several assemblies together.
    Acts as a new assembly.
    Newbler is outdated as Roche has stopped making sequencing machines in 2013.
    Newbler is closed source and Roche dismissed the community's
    petition to make it open source.
    You can only obtain the software via a registration form."""

    short_name = 'newbler'
    long_name  = 'Newbler from Roche GS De Novo Assembler software v2.9'
    executable = 'runAssembly'
    url        = ''
    dependencies = []

    all_paths = """
    /combined_cut_up.fasta
    /output/454AllContigs.fna
    /stdout.txt
    /stderr.txt
    /filtered_contigs.fasta
    /bins/
    /report/report.pdf
    /hit_profile/
    /traits/
    /bins_summary/
    """

    def __repr__(self): return '<%s object kmer %s>' % (self.__class__.__name__, self.kmer_size)
    def __nonzero__(self): return bool(self.p.filtered_contigs)
    def __len__(self):  return len(self.samples)
    def __getitem__(self, key):
        if isinstance(key, basestring): return [c for c in self.children if c.name == key][0]
        return self.children[key]

    def __init__(self, samples, assemblies, result_dir, slice_length=1000):
        # Base parameters #
        self.samples      = samples
        self.assemblies   = assemblies
        self.children     = samples
        self.result_dir   = result_dir
        self.slice_length = slice_length
        # Auto paths #
        self.base_dir = self.result_dir + self.short_name + '/'
        self.p = AutoPaths(self.base_dir, self.all_paths)
        # Report #
        self.report = AssemblyReport(self)
        # Name #
        self.name = self.short_name + '_' + str(len(self.assemblies))

    def run(self, cpus=None, verbose=True):
        # Variable threads #
        if cpus is None: cpus = num_processors
        # Check that all the sets of contigs have a cut-up version #
        for assembly in self.assemblies:
            if not assembly.p.cut_up:
                if verbose: print 'Cutting up contigs from %s' % assembly; sys.stdout.flush()
                params = (assembly.results.contigs_fasta, assembly.p.cut_up)
                shell_output("fasta_cut_up %s > %s" % params)
        # Merge the cut-uped fastas #
        if not self.p.cut_up:
            if verbose: print 'Combining contigs from %s' % self.assemblies; sys.stdout.flush()
            paths = [assembly.p.cut_up.path for assembly in self.assemblies]
            shell_output('cat %s > %s' % (' '.join(paths), self.p.combined))
        # Call newbler #
        if verbose: print 'Calling Newbler...'; sys.stdout.flush()
        sh.runAssembly('-force',
                       '-cpu', cpus,
                       '-o', self.p.output_dir.path,
                       self.p.combined.path,
                       _out=self.p.stdout.path,
                       _err=self.p.stderr.path)
        # Check there is something #
        contigs = FASTA(self.p.fna)
        if len(contigs) == 0: raise Exception("Newbler found exactly 0 contigs in your dataset.")
        # Filter short contigs #
        filtered = FASTA(self.p.filtered)
        contigs.extract_length(new_path=filtered, lower_bound=self.slice_length)
        # Make indexes (used later) #
        if not os.path.exists(filtered + '.1.bt2'): filtered.index_bowtie()
        if not os.path.exists(filtered + '.fai'):   filtered.index_samtools()

    @property
    def short_description(self):
        return "Newbler"

    @property
    def description(self):
        return "Newbler merging %i assemblies (%s)" \
            % (len(self.assemblies), self.kmer_size)

    @property
    def kmer_size(self): return andify(str(a.kmer_size) for a in self.assemblies)
    @property
    def length_cutoff(self): return andify(list(set(str(a.length_cutoff) for a in self.assemblies)))

    @property_cached
    def results(self):
        results = NewblerResults(self)
        if not results: raise Exception("You can't access results from Newbler before running the algorithm.")
        return results

###############################################################################
class NewblerResults(AssemblyResults):

    @property_cached
    def mappings(self):
        """Map each of the samples used in the assembly back to this assembly."""
        return OrderedDict([(s.name, s.mapper_merged) for s in self.newbler.samples])
