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
mem_size = nr_threads * 3

###############################################################################
class Mapper(object):
    """Maps reads from a Pool to an Assembly."""

    all_paths = """
    /map.sam
    /map.bam
    /map_s.bam
    /map_smd.bam
    /map_smds.bam
    /map_smd.metrics
    /map_smds.coverage
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
        """Maps reads from self.pool to self.assembly using bowtie2. PCR
        Duplicates are afterwards removed using MarkDuplicates. BEDTools is
        used to determine coverage."""
        # Create bowtie2 assembly index #
        contigs = self.assembly.contigs_fasta
        sh.bowtie2_build(contigs, contigs)
        # Create the mapping #
        sh.bowtie2('-p', nr_threads, '-x', contigs, '-1', self.pool.fwd, '-2', self.pool.rev, '-S', self.p.sam)
        # Create samtools assembly index, sort and index bamfile #
        sh.samtools('faidx', contigs)
        sh.samtools('view', '-bt', contigs + '.fai', self.p.sam, _out=self.p.map_bam)
        sh.samtools('sort', self.p.map_bam, self.p.map_s_bam)
        sh.samtools('index', self.p.map_s_bam)
        # Remove PCR duplicates with MarkDuplicates #
        sh.java('-Xms1g', '-Xmx' + str(mem_size) + 'g', '-XX:ParallelGCThreads=' +
                str(nr_threads), '-XX:MaxPermSize=1g',
                '-XX:+CMSClassUnloadingEnabled', '-jar', 'MarkDuplicates.jar',
                'INPUT=' + self.p.map_s_bam, 'OUTPUT=' + self.p.map_smd_bam,
                'METRICS=' + self.p.map_smds_metrics, 'AS=TRUE',
                'VALIDATION_STRINGENCY=LENIENT',
                'MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000',
                'REMOVE_DUPLICATES=TRUE')
        # Sort and index bam without duplicates #
        sh.samtools('sort', self.p.map_smd_bam, self.p.map_smds_bam)
        sh.samtools('index', self.p.map_smds_bam)
        # Determine Coverage with BEDTools #
        sh.genomeCoverageBed('-ibam', self.p.map_smds_bam, _out=self.p.map_smds_coverage)
