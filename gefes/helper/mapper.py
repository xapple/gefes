# Futures #
from __future__ import division

# Built-in modules #
import os

# Internal modules #
from gefes.common.autopaths import AutoPaths

# Third party modules #
import sh

# Constant #
nr_threads = os.environ['SLURM_JOB_CPUS_PER_NODE']

###############################################################################
class Mapper(object):
    """Lorem ipsum."""

    all_paths = """
    /map.sam
    /map.bam
    """

    def __repr__(self): return '<%s object of %s>' % (self.__class__.__name__, self.parent)

    def __init__(self, pool, assembly):
        # Save attributes #
        self.parent, self.pool = pool, pool
        self.assembly = assembly
        # Auto paths #
        self.base_dir = self.pool.p.mapping_dir + self.assembly.short_name
        self.p = AutoPaths(self.base_dir, self.all_paths)

    def map(self):
        # Create bowtie2 assembly index #
        contigs = self.assembly.contigs_fasta
        sh.bowtie2_build(contigs, contigs)
        # Create the mapping #
        sh.bowtie2('-p', nr_threads, '-x', contigs, '-1', self.pool.fwd, '-2', self.pool.rev, '-S', self.p.sam)
        # Create samtools assembly index #
        sh.samtools('faidx', self.assembly.contigs_fasta.path)
        sh.samtools('view', '-bt', self.assembly.contigs_fasta.path + '.fai', self.p.sam, _out=self.p.bam)
