# Futures #
from __future__ import division

# Built-in modules #
import os

# Internal modules #
from gefes.common.autopaths import AutoPaths

# Third party modules #
import sh

###############################################################################
class Mapper(object):
    """Lorem ipsum."""

    all_paths = """
    /map/map.sam
    /map/map.bam
    """

    def __repr__(self): return '<%s object of %s>' % (self.__class__.__name__, self.parent)

    def __init__(self, pool, assembly):
        # Save parents #
        self.parent, self.pool = pool, pool
        self.assembly = assembly
        # Auto paths #
        self.base_dir = self.assembly.base_dir + self.pool.id_name + '/'
        self.p = AutoPaths(self.base_dir, self.all_paths)
        # Convenience objects #
        self.fwd = self.pool.fwd
        self.rev = self.pool.rev

    def map(self):
        # Create bowtie2 assembly index #
        ref = self.assembly.contigs_fasta.path

        if not os.path.exists(ref + ".1.bt2"):
            sh.bowtie2_build(ref, ref)
        nr_threads = 1
        sh.bowtie2('-p', nr_threads, '-x', ref, '-1', self.fwd.path, '-2', self.rev.path, '-S', self.p.sam)

        # Create samtools assembly index #
        if not os.path.exists(self.assembly.contigs_fasta.path + ".fai"):
            sh.samtools('faidx', self.assembly.contigs_fasta.path)

        sh.samtools('view', '-bt', self.assembly.contigs_fasta.path + '.fai', self.p.sam, _out=self.p.bam)
