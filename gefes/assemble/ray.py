# Futures #
from __future__ import division

# Built-in modules #
import os, socket

# Internal modules #
from gefes.assemble.contig import Contig
from plumbing.common import flatter
from plumbing.autopaths import AutoPaths
from plumbing.cache import property_cached
from plumbing.slurm import nr_threads
from fasta import FASTA

# Third party modules #
import sh

# Constant #
hostname = socket.gethostname()

###############################################################################
class Ray(object):
    """Will run the co-assembly of several samples by calling the Ray assembler.
    https://github.com/sebhtml/ray"""

    short_name = 'ray'
    executable = 'ray231'

    all_paths = """
    /stdout.txt
    /stderr.txt
    /output/
    /output/Contigs.fasta
    /output/report.txt
    """

    def __repr__(self): return '<%s object of %s>' % (self.__class__.__name__, self.parent)

    def __init__(self, samples, result_dir, kmer_size=41):
        # Base params #
        self.samples = samples
        self.result_dir = result_dir
        self.kmer_size = kmer_size
        self.base_dir = self.result_dir + 'ray/%i/' % self.kmer_size
        self.p = AutoPaths(self.base_dir, self.all_paths)
        # To be set when run #
        self.out_dir = None
        self.paths = []

    @property_cached
    def contigs_fasta(self): return FASTA(self.p.Contigs)

    @property_cached
    def contigs(self):
        """A list of all the contigs produced as custom objects"""
        return [Contig(self, s) for s in self.contigs_fasta]

    def run(self):
        # Check samples #
        for s in self.samples:
            if not s.loaded: s.load()
            assert s.clean.exists
            assert s.singletons
        # Ray needs a non-existing directory #
        self.out_dir = self.p.output_dir
        self.out_dir.remove()
        # Make the pairs of fastq #
        self.paths = lambda s: ('-p', s.clean.fwd.path, s.clean.rev.path, '-s', s.singletons.path)
        self.paths = flatter([self.paths(s) for s in self.samples])
        # Call Ray on different setting #
        if os.environ.get('CSCSERVICE') == 'sisu': stats = self.sisu()
        elif os.environ.get('SLURM_JOB_PARTITION') == 'halvan': stats = self.halvan()
        elif os.environ.get('SNIC_RESOURCE') == 'milou': stats = self.milou()
        else: stats = self.local()
        # Print the report #
        with open(self.p.report, 'w') as handle: handle.write(str(stats))

    def sisu(self): return sh.aprun('-n', nr_threads, self.executable, '-k', self.kmer_size, '-o', self.out_dir, *self.paths, _out=self.p.stdout.path, _err=self.p.stderr.path)

    def halvan(self): return sh.mpiexec('-n', nr_threads, self.executable, '-k', self.kmer_size, '-o', self.out_dir, *self.paths, _out=self.p.stdout.path, _err=self.p.stderr.path)

    def milou(self):
        if 'SLURM_NODELIST' not in os.environ: os.environ['SLURM_NODELIST'] = hostname
        commands = ['-n', nr_threads, self.executable, '-k', self.kmer_size, '-o', self.out_dir] + self.paths
        return sh.mpiexec(*commands, _out=self.p.stdout.path, _err=self.p.stderr.path)

    def local(self):
        ray = sh.Command(self.executable)
        return ray('-k', self.kmer_size, '-o', self.out_dir, *self.paths, _out=self.p.stdout.path, _err=self.p.stderr.path)

###############################################################################
class RayResults(object):

    all_paths = """
    /lorem
    """
    def __nonzero__(self): return self.dest.exists

    def __init__(self, dest):
        self.dest = dest