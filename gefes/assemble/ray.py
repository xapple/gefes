# Futures #
from __future__ import division

# Built-in modules #
import os, socket

# Internal modules #
from gefes.helper.contig import Contig
from plumbing.common import flatten
from plumbing.autopaths import AutoPaths
from plumbing.cache import property_cached
from fasta import FASTA
from plumbing.slurm import nr_threads

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
    kmer_size = 51

    all_paths = """
    /graphs/
    /ray_output/frame.csv
    /ray_output/
    /ray_output/Contigs.fasta
    /ray_output/report.txt
    /metapathways/
    """

    def __repr__(self): return '<%s object of %s>' % (self.__class__.__name__, self.parent)

    def __init__(self, aggregate):
        # Save parent #
        self.parent, self.aggregate = aggregate, aggregate
        # Auto paths #
        self.base_dir = self.parent.p.assembly_dir
        self.p = AutoPaths(self.base_dir, self.all_paths)
        # Convenience objects #
        self.contigs_fasta = FASTA(self.p.Contigs)

    @property_cached
    def contigs(self):
        """A list of all the contigs produced as custom objects"""
        return [Contig(self, s) for s in self.contigs_fasta]

    def run(self):
        # Ray needs a non-existing directory #
        out_dir = self.p.output_dir
        out_dir.remove()
        # Make the pairs of fastq #
        pairs = flatten([('-p', p.cleaner.fwd.path, p.cleaner.rev.path) for p in self.parent])
        # Call Ray on the cray #
        if os.environ.get('CSCSERVICE') == 'sisu':
            stats = sh.aprun('-n', nr_threads, self.executable, '-k', self.kmer_size, '-o', out_dir, *pairs)
        # Call Ray on Kalkyl #
        elif os.environ.get('SNIC_RESOURCE') == 'kalkyl':
            stats = sh.mpiexec('-n', nr_threads, self.executable, '-k', self.kmer_size, '-o', out_dir, *pairs)
        # Call Ray just locally #
        else:
            ray = sh.Command(self.executable)
            stats = ray('-k', self.kmer_size, '-o', out_dir, *pairs)
        # Print the report #
        with open(self.p.report, 'w') as handle: handle.write(str(stats))

###############################################################################
class RayResults(object):

    all_paths = """
    /lorem
    """
    def __nonzero__(self): return self.dest.exists

    def __init__(self, dest):
        self.dest = dest