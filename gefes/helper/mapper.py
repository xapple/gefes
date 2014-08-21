# Futures #
from __future__ import division

# Built-in modules #
import os

# Internal modules #
import gefes
from gefes.common.autopaths import AutoPaths
from gefes.common.cache import property_cached
from gefes.helper.linkage import parse_linkage_info_bam
from gefes.common.slurm import nr_threads

# Third party modules #
import sh, pandas, os

###############################################################################
class Mapper(object):
    """Maps reads from a Pool object to an Assembly object.
    Names follow this standard:
      * The 'map_s' file is sorted.
      * The 'map_smd' file is sorted and mark duplicated.
      * The 'map_smd' file is sorted, duplicated, and sorted again."""

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

    def map(self,threads=nr_threads):
        """Maps reads from self.pool to self.assembly using bowtie2. PCR
        Duplicates are afterwards removed using MarkDuplicates. BEDTools is
        used to determine coverage."""
        # Check indexes #
        if not os.path.exists(self.contigs + '.1.bt2'):
            raise(Exception('Bowtie2 index file not created, run index() first'))
        if not os.path.exists(self.contigs + '.fai'):
            raise(Exception('Samtools index file not created, run index() first'))
        # Do the mapping #
        sh.bowtie2('-p', threads, '-x', self.contigs, '-1', self.pool.fwd, '-2', self.pool.rev, '-S', self.p.map_sam)
        # Create bam, sort and index bamfile #
        sh.samtools('view', '-bt', self.contigs + '.fai', self.p.map_sam, '-o', self.p.map_bam)
        sh.samtools('sort', self.p.map_bam, self.p.map_s_bam.prefix_path)
        sh.samtools('index', self.p.map_s_bam)
        # Remove PCR duplicates #
        self.remove_duplicates()
        # Sort and index bam without duplicates #
        sh.samtools('sort', self.p.map_smd_bam, self.p.map_smds_bam.prefix_path)
        sh.samtools('index', self.p.map_smds_bam)
        # Compute coverage #
        sh.genomeCoverageBed('-ibam', self.p.map_smds_bam, _out=str(self.p.map_smds_coverage)) 
        # Clean up #
        os.remove(self.p.map_sam)
        os.remove(self.p.map_bam)
        os.remove(self.p.map_smd_bam)

    def remove_duplicates(self):
        """Remove PCR duplicates with MarkDuplicates."""
        # Estimate size #
        mem_size = "4"
        
        sh.java('-Xmx%sg' % mem_size,
                '-XX:ParallelGCThreads=%s' % nr_threads,
                '-XX:+CMSClassUnloadingEnabled',
                '-jar', gefes.repos_dir + 'bin/picard-tools-1.101/MarkDuplicates.jar',
                'INPUT=%s' % self.p.map_s_bam,
                'OUTPUT=%s' % self.p.map_smd_bam,
                'METRICS_FILE=%s' % self.p.map_smd_metrics,
                'AS=TRUE',
                'VALIDATION_STRINGENCY=LENIENT',
                'MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000',
                'REMOVE_DUPLICATES=TRUE')

    #-------------------------------------------------------------------------#
    @property_cached
    def coverage(self):
        """Uses the BEDTools genomeCoverageBed histogram output to determine mean
        coverage and percentage covered for each contig. Returns a dict with fasta id as
        key and containing 'percent_covered' with 'cov_mean' information inside.
        The output file has the following headers:
          * headers = ['name', 'depth', 'count', 'length', 'fraction']
        http://bedtools.readthedocs.org/en/latest/content/tools/genomecov.html"""
        # Main loop #
        out_dict = {}
        with open(self.p.map_smds_coverage) as handler:
            for line in handler:
                name, depth, count, length, fraction = line.split()
                d = out_dict.setdefault(name, {"cov_mean": 0, "percent_covered": 100})
                if int(depth) == 0:
                    d["percent_covered"] = 100 - float(fraction) * 100.0
                else:
                    d["cov_mean"] = d.get("cov_mean", 0) + int(depth) * float(fraction)
        # Add 0 coverage for contigs not in the output #
        for c in self.assembly.contigs:
            if c.name not in out_dict:
                out_dict[c.name] = {"cov_mean": 0, "percent_covered": 0}
        # Return result #
        return out_dict

    @property_cached
    def coverage_alternative(self):
        """Another way of computing the same thing. Not functional yet."""
        frame = pandas.io.parsers.read_csv(self.p.map_smds_coverage, sep='\t')
        out_dict = {}
        for name in self.assembly.contigs.names:
            cov_mean = frame[name][3] * frame[name][4] if name in frame else 0
            percent_covered = 100 if 0 not in frame[name][3] else 100 - frame[name][3][0] * 100.0
            out_dict[name] = {"cov_mean": cov_mean, "percent_covered": percent_covered}
        return out_dict

    #-------------------------------------------------------------------------#
    @property_cached
    def linkage_and_readcount(self):
        return parse_linkage_info_bam(bamfile=self.p.map_smds_bam.path,
                readlength=100, min_contig_length=100, regionlength=500,
                fullsearch=False)

    @property
    def linkage(self):
        return self.linkage_and_readcount[0]

    @property
    def readcount(self):
        return self.linkage_and_readcount[1]
