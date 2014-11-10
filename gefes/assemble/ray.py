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
    """

    def __repr__(self): return '<%s object kmer %i>' % (self.__class__.__name__, self.kmer_size)

    def __init__(self, samples, result_dir, kmer_size=41):
        # Base parameters #
        self.samples = samples
        self.result_dir = result_dir
        self.kmer_size = kmer_size
        self.base_dir = self.result_dir + 'ray/%i/' % self.kmer_size
        self.p = AutoPaths(self.base_dir, self.all_paths)
        # To be set when run #
        self.out_dir = None
        self.paths = []

    def run(self):
        # Check samples #
        for s in self.samples:
            if not s.loaded: s.load()
            assert s.clean.exists
            assert s.singletons
        # Ray needs a non-existing directory otherwise it is unhappy #
        self.out_dir = self.p.output_dir
        self.out_dir.remove()
        # Make the pairs of fastq #
        self.paths = lambda s: ('-p', s.clean.fwd.path, s.clean.rev.path, '-s', s.singletons.path)
        self.paths = flatter([self.paths(s) for s in self.samples])
        # Call Ray on different setting #
        if os.environ.get('CSCSERVICE') == 'sisu':              stats = self.sisu()
        elif os.environ.get('SLURM_JOB_PARTITION') == 'halvan': stats = self.halvan()
        elif os.environ.get('SNIC_RESOURCE') == 'milou':        stats = self.milou()
        else:                                                   stats = self.local()
        # Print the report #
        with open(self.p.report, 'w') as handle: handle.write(str(stats))

    def sisu(self):
        """Run the assembly on the Sisu super computer at sisu-login1.csc.fi"""
        return sh.aprun('-n', nr_threads, self.executable, '-k', self.kmer_size, '-o', self.out_dir, *self.paths, _out=self.p.stdout.path, _err=self.p.stderr.path)

    def halvan(self):
        """Run the assembly on the large memory computer at http://www.uppmax.uu.se/halvan-user-guide"""
        return sh.mpiexec('-n', nr_threads, self.executable, '-k', self.kmer_size, '-o', self.out_dir, *self.paths, _out=self.p.stdout.path, _err=self.p.stderr.path)

    def milou(self):
        """Run the assembly on one node of the milou cluster at http://www.uppmax.uu.se/milou-user-guide"""
        if 'SLURM_NODELIST' not in os.environ: os.environ['SLURM_NODELIST'] = hostname
        commands = ['-n', nr_threads, self.executable, '-k', self.kmer_size, '-o', self.out_dir] + self.paths
        return sh.mpiexec(*commands, _out=self.p.stdout.path, _err=self.p.stderr.path)

    def local(self):
        """Run the assembly on the computer where this pipeline is running"""
        ray = sh.Command(self.executable)
        return ray('-k', self.kmer_size, '-o', self.out_dir, *self.paths, _out=self.p.stdout.path, _err=self.p.stderr.path)

    @property_cached
    def results(self):
        results = RayResults(self)
        if not results: raise Exception("You can't access results from Ray before running the assembly.")
        return results

###############################################################################
class RayResults(object):

    all_paths = """
    /output/Contigs.fasta
    /output/report.txt
    """

    def __nonzero__(self): return self.contigs_fasta.exists
    def __init__(self, ray):
        self.ray = ray
        self.p = AutoPaths(self.ray.base_dir, self.all_paths)

    @property_cached
    def contigs_fasta(self):
        """Just the contigs as a FASTA file object"""
        return FASTA(self.p.Contigs)

    @property_cached
    def contigs(self):
        """A list of all the contigs produced returned as our Contig custom objects"""
        return [Contig(self.ray, s) for s in self.contigs_fasta]
