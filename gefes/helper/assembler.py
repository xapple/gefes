# Futures #
from __future__ import division

# Built-in modules #
import os

# Internal modules #
from gefes.common import flatten
from gefes.common.autopaths import AutoPaths
from gefes.common.cache import property_cached
from gefes.graphs import assembly_plots
from gefes.helper.contig import Contig
from gefes.fasta.single import FASTA

# Third party modules #
import sh

# Constant #
nr_threads = os.environ['SLURM_JOB_CPUS_PER_NODE']

###############################################################################
class Assembly(object):
    """The co-assembly of all pools."""
    short_name = 'ray'

    all_paths = """
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

    @property_cached
    def contigs(self):
        return [Contig(self, s) for s in self.contigs_fasta]

    def assemble(self):
        # Ray needs a non-existing directory #
        self.p.output_dir.remove()
        # Call Ray #
        stats = sh.mpiexec('-n', nr_threads, 'Ray', '-k', 81, '-o', self.p.output_dir,
            *flatten([('-p', p.cleaner.fwd.path, p.cleaner.rev.path) for p in self.parent]))
        # Print the report #
        with open(self.p.report, 'w') as handle: handle.write(str(stats))

    def make_plots(self):
        for cls_name in assembly_plots.__all__:
            cls = getattr(assembly_plots, cls_name)
            cls(self).plot()
