# Futures #
from __future__ import division

# Built-in modules #
import os

# Internal modules #
from gefes.common.autopaths import AutoPaths
from gefes.common.cache import property_cached
from gefes.helper.linkage import parse_linkage_info_bam
import gefes

# Third party modules #
import sh

# Constant #
nr_threads = int(os.environ.get('SLURM_JOB_CPUS_PER_NODE', 1))

###############################################################################
class Mapper(object):
    """Maps reads from a Pool to an Assembly.
    map_s is sorted.
    map_smd is sorted and mark duplicated."""

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
        self.contigs = self.assembly.contigs_fasta
        # Auto paths #
        self.base_dir = self.pool.p.mapping_dir + self.assembly.short_name
        self.p = AutoPaths(self.base_dir, self.all_paths)

    def map(self):
        """Maps reads from self.pool to self.assembly using bowtie2. PCR
        Duplicates are afterwards removed using MarkDuplicates. BEDTools is
        used to determine coverage."""
        # Check indexes #
        if not os.path.exists(self.contigs + '.1.bt2'):
            raise(Exception('Bowtie2 index file not created, run index_assembly()'))
        if not os.path.exists(self.contigs + '.fai'):
            raise(Exception('Samtools index file not created, run index_assembly()'))
        # Do the mapping #
        sh.bowtie2('-p', nr_threads, '-x', self.contigs, '-1', self.pool.fwd, '-2', self.pool.rev, '-S', self.p.sam)
        # Create bam, Sort and index bamfile #
        sh.samtools('view', '-bt', self.contigs + '.fai', self.p.sam, _out=self.p.map_bam)
        sh.samtools('sort', self.p.map_bam, self.p.map_s_bam[:-4])
        sh.samtools('index', self.p.map_s_bam)
        # Remove PCR duplicates #
        self.remove_duplicates()
        # Sort and index bam without duplicates #
        sh.samtools('sort', self.p.map_smd_bam, self.p.map_smds_bam[:-4])
        sh.samtools('index', self.p.map_smds_bam)
        self.calc_coverage()
        # Clean up #
        os.remove(self.p.sam)
        os.remove(self.p.map_bam)
        os.remove(self.p.map_smd_bam)

    def remove_duplicates(self):
        """Remove PCR duplicates with MarkDuplicates."""
        # Estimate size #
        mem_size = nr_threads * 2
        perm_size = str(max(int(mem_size * 0.3), 1))
        sh.java('-Xms%sg' % perm_size,
                '-Xmx%sg' % mem_size,
                '-XX:ParallelGCThreads=%s' % nr_threads,
                '-XX:MaxPermSize=%sg' % perm_size,
                '-XX:+CMSClassUnloadingEnabled',
                '-jar', gefes.repos_dir + 'bin/picard-tools-1.101/MarkDuplicates.jar',
                'INPUT=%s' % self.p.map_s_bam,
                'OUTPUT=%s' % self.p.map_smd_bam,
                'METRICS_FILE=%s' % self.p.map_smd_metrics,
                'AS=TRUE',
                'VALIDATION_STRINGENCY=LENIENT',
                'MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000',
                'REMOVE_DUPLICATES=TRUE')

    def calc_coverage(self):
        # Determine Coverage with BEDTools #
        sh.genomeCoverageBed('-ibam', self.p.map_smds_bam, _out=self.p.map_smds_coverage)

    @property_cached
    def coverage(self):
        """Uses the BEDTools genomeCoverageBed histogram output to determine mean
        coverage and percentage covered for each contig.  Returns dict with fasta id as
        key and percentage covered and cov_mean as keys for the inner dictionary."""
        out_dict = {}

        with open(self.p.map_smds_coverage) as handler:
            for line in handler:
                name, depth, count, length, fraction = line.split()

                d = out_dict.setdefault(name, {"cov_mean": 0, "percentage_covered": 100})

                if int(depth) == 0:
                    d["percentage_covered"] = 100 - float(fraction) * 100.0
                else:
                    d["cov_mean"] = d.get("cov_mean", 0) + int(depth) * float(fraction)

        return out_dict

    @property_cached
    def linkage_and_readcount(self):
        return parse_linkage_info_bam(bamfile=self.p.map_smds_bam.path,
                readlength=100, min_contig_length=100, regionlength=500,
                fullsearch=True)

    @property
    def linkage(self):
        return self.linkage_and_readcount[0]

    @property
    def readcount(self):
        return self.linkage_and_readcount[1]
