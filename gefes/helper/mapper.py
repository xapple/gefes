# Futures #
from __future__ import division

# Built-in modules #
import os

# Internal modules #
from gefes.common.autopaths import AutoPaths
import gefes

# Third party modules #
import sh

# Constant #
try:
    nr_threads = int(os.environ['SLURM_JOB_CPUS_PER_NODE'])
except KeyError:
    nr_threads = 1
mem_size = nr_threads * 2

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

    def index_assembly(self):
        """Create index for both bowtie2 and samtools on assembly fasta file."""
        contigs = self.assembly.contigs_fasta
        sh.bowtie2_build(contigs, contigs)
        sh.samtools('faidx', contigs)

    def remove_duplicates(self):
        """Remove PCR duplicates with MarkDuplicates."""
        sh.java('-Xms' + str(max(int(mem_size * 0.3), 1)) + 'g', '-Xmx' +
                str(mem_size) + 'g', '-XX:ParallelGCThreads=' +
                str(nr_threads), '-XX:MaxPermSize=' + str(max(int(mem_size * 0.3), 1)) + 'g',
                '-XX:+CMSClassUnloadingEnabled', '-jar', gefes.repos_dir +
                'bin/picard-tools-1.101/MarkDuplicates.jar', 'INPUT=' +
                self.p.map_s_bam, 'OUTPUT=' + self.p.map_smd_bam,
                'METRICS_FILE=' + self.p.map_smd_metrics, 'AS=TRUE',
                'VALIDATION_STRINGENCY=LENIENT',
                'MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000',
                'REMOVE_DUPLICATES=TRUE')

    def calc_coverage(self):
        # Determine Coverage with BEDTools #
        sh.genomeCoverageBed('-ibam', self.p.map_smds_bam, _out=self.p.map_smds_coverage)

    def map(self):
        """Maps reads from self.pool to self.assembly using bowtie2. PCR
        Duplicates are afterwards removed using MarkDuplicates. BEDTools is
        used to determine coverage."""
        # Create bowtie2 assembly index #
        contigs = self.assembly.contigs_fasta
        if not os.path.exists(contigs + '.1.bt2'):
            raise(Exception('Bowtie2 index file not created, run index_assembly()'))
        # Do the mapping #
        sh.bowtie2('-p', nr_threads, '-x', contigs, '-1', self.pool.fwd, '-2', self.pool.rev, '-S', self.p.sam)
        # Create bam, Sort and index bamfile #
        if not os.path.exists(contigs + '.fai'):
            raise(Exception('Samtools index file not created, run index_assembly()'))
        sh.samtools('view', '-bt', contigs + '.fai', self.p.sam, _out=self.p.map_bam)
        sh.samtools('sort', self.p.map_bam, self.p.map_s_bam[:-4])
        sh.samtools('index', self.p.map_s_bam)
        # Remove PCR duplicates #
        self.remove_duplicates()
        # Sort and index bam without duplicates #
        sh.samtools('sort', self.p.map_smd_bam, self.p.map_smds_bam[:-4])
        sh.samtools('index', self.p.map_smds_bam)
        self.calc_coverage()
