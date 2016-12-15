# Futures #
from __future__ import division

# Built-in modules #
import os, warnings

# Internal modules #
from gefes.assemble        import Assembler, AssemblyResults
from gefes.report.assembly import AssemblyReport

# First party modules #
from plumbing.autopaths import AutoPaths
from plumbing.cache import property_cached
from plumbing.slurm import num_processors, current_server
from fasta import FASTA

# Third party modules #
import sh

###############################################################################
class Megahit(Assembler):
    """Will run the co-assembly of several samples by calling the Megahit assembler.
    Expects version 1.0.6.1
    We remove all the contigs below the length cutoff threshold."""

    short_name = 'megahit'
    long_name  = 'MEGAHIT assembler 1.0.6.1'
    executable = 'megahit'
    url        = 'https://github.com/voutcn/megahit'
    dependencies = []

    all_paths = Assembler.all_paths + """
    /output/final.contigs.fa
    /output/log
    /output/opts.txt
    /filtered.fasta
    /stdout.txt
    /stderr.txt
    """

    kmer_size = 'variable'

    def __init__(self, samples, result_dir, length_cutoff=1000):
        # Base parameters #
        self.samples       = samples
        self.children      = samples
        self.result_dir    = result_dir
        self.length_cutoff = length_cutoff
        # Auto paths #
        self.base_dir = self.result_dir + self.short_name + '/'
        self.p = AutoPaths(self.base_dir, self.all_paths)

    @property_cached
    def report(self):
        """The PDF report."""
        return AssemblyReport(self)

    def run(self, cpus=None):
        # Variable threads #
        if cpus is None: cpus = num_processors
        # Check samples #
        for s in self.samples:
            assert s.clean.exists
            assert s.singletons
        # Remove outdir #
        self.out_dir = self.p.output_dir
        self.out_dir.remove()
        # Make the pairs of fastq #
        self.paths  = ['-1', ','.join(s.clean.fwd.path  for s in self.samples)]
        self.paths += ['-2', ','.join(s.clean.rev.path  for s in self.samples)]
        self.paths += ['-r', ','.join(s.singletons.path for s in self.samples)]
        # Executable #
        megahit = sh.Command(self.executable)
        # Check version #
        assert "v1.0.6.1" in megahit('--version')
        # Call it #
        megahit('--out-dir',           self.out_dir,
                '--presets',          'meta-large',
                '--memory',           '0.7',
                '--num-cpu-threads',   cpus,
                '--tmp-dir',          '/tmpfs/lucas/tmp/megahit/',
                '--verbose',
                *self.paths,
                _out=self.p.stdout.path,
                _err=self.p.stderr.path)
        # Check it worked #
        if not self.p.contigs.exists: raise Exception("Did not create any contigs.")
        # Check there is something #
        contigs = FASTA(self.p.contigs)
        if len(contigs) == 0: raise Exception("Found exactly 0 contigs in your dataset.")
        # Filter short contigs #
        filtered = FASTA(self.p.filtered)
        contigs.extract_length(new_path=filtered, lower_bound=self.length_cutoff)
        message = "After filtering, there were no contigs left in sample '%s'" % self
        if not filtered: warnings.warn(message)
        # Make indexes (used later) #
        if filtered and not os.path.exists(filtered + '.1.bt2'): filtered.index_bowtie()
        if filtered and not os.path.exists(filtered + '.fai'):   filtered.index_samtools()

    @property
    def short_description(self):
        return "Megahit with kmer %i " % self.kmer_size

    @property
    def description(self):
        return "Megahit with kmer %i and %i bp cutoff" % (self.kmer_size, self.length_cutoff)

    @property_cached
    def results(self):
        results = MegahitResults(self)
        if not results: raise Exception("You can't access results from Megahit before running the assembly.")
        return results

###############################################################################
class MegahitResults(AssemblyResults):
    pass