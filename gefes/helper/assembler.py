# Futures #
from __future__ import division

# Built-in modules #
import os, socket

# Internal modules #
from gefes.common import flatten
from gefes.common.autopaths import AutoPaths
from gefes.common.cache import property_cached
from gefes.graphs import assembly_plots
from gefes.helper.contig import Contig
from gefes.fasta.single import FASTA
from gefes.common.slurm import nr_threads

# Third party modules #
import sh

# Constant #
hostname = socket.gethostname()

###############################################################################
class Assembly(object):
    """The co-assembly of all pools."""
    short_name = 'ray'

    all_paths = """
    /graphs/
    /ray_output/
    /ray_output/Contigs.fasta
    /ray_output/report.txt
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
        # Graphs #
        self.graphs = [getattr(assembly_plots, cls_name)(self) for cls_name in assembly_plots.__all__]

    @property_cached
    def contigs(self):
        return [Contig(self, s) for s in self.contigs_fasta]

    def assemble(self):
        # Ray needs a non-existing directory #
        out_dir = self.p.output_dir
        out_dir.remove()
        # Make the pairs of fastq #
        pairs = flatten([('-p', p.cleaner.fwd.path, p.cleaner.rev.path) for p in self.parent])
        # Call Ray on the cray #
        if os.environ.get('CSCSERVICE') == 'sisu':
            stats = sh.aprun('-n', nr_threads, 'Ray23', '-k', 81, '-o', out_dir, *pairs)
        # Call Ray on Kalkyl #
        elif os.environ.get('SNIC_RESOURCE') == 'kalkyl':
            stats = sh.mpiexec('-n', nr_threads, 'Ray23', '-k', 81, '-o', out_dir, *pairs)
        # Call Ray just locally #
        else:
            stats = sh.Ray23('-k', 81, '-o', out_dir, *pairs)
        # Print the report #
        with open(self.p.report, 'w') as handle: handle.write(str(stats))

    def index(self):
        """Create two indexes. For both bowtie2 and samtools on assembly fasta file."""
        sh.bowtie2_build(self.contigs_fasta, self.contigs_fasta)
        sh.samtools('faidx', self.contigs_fasta)

    def make_plots(self):
        for graph in self.graphs: graph.plot()