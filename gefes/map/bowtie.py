# Futures #
from __future__ import division

# Built-in modules #
import os

# Internal modules #
import gefes

# First party modules #
from plumbing.autopaths import AutoPaths
from plumbing.cache import property_cached, property_pickled
from plumbing.slurm import num_processors

# Third party modules #
import sh, pandas

###############################################################################
class Bowtie(object):
    """Use Bowtie2 at http://bowtie-bio.sourceforge.net/bowtie2/index.shtml
    to maps reads from a Sample object back to an Assembly object.
    Expects version 2.2.4.
    SAMtools is used to index and sort the result (v0.1.19).
    PCR Duplicates are subsequently removed using MarkDuplicates (v1.101).
    BEDTools is then used to determine coverage (v2.15.0).

    Names follow this standard:
      * The 'map_s' file is sorted.
      * The 'map_smd' file is sorted and mark duplicated.
      * The 'map_smd' file is sorted, duplicated, and sorted again.
    """

    short_name = 'bowtie'

    all_paths = """
    /map.sam
    /map.bam
    /map_s.bam
    /map_smd.bam
    /map_smds.bam
    /map_smds.bai
    /map_smd.metrics
    /map_smds.coverage
    /statistics.pickle
    """

    def __repr__(self): return '<%s object of %s on %s>' % \
                        (self.__class__.__name__, self.sample, self.assembly)

    def __init__(self, sample, assembly, result_dir):
        # Save attributes #
        self.sample = sample
        self.assembly = assembly
        self.result_dir = result_dir
        # Convenience shortcuts #
        self.contigs_fasta = self.assembly.results.contigs_fasta
        # Auto paths #
        self.base_dir = self.result_dir + self.short_name + '/'
        self.p = AutoPaths(self.base_dir, self.all_paths)

    def run(self):
        # Check both type of indexes exist #
        if not os.path.exists(self.contigs_fasta + '.1.bt2'): self.contigs_fasta.index_bowtie()
        if not os.path.exists(self.contigs_fasta + '.fai'): self.contigs_fasta.index_samtools()
        # Make our options #
        options = ['-p', num_processors,
                   '-x', self.assembly.results.contigs_fasta,
                   '-1', self.sample.fwd_path,
                   '-2', self.sample.rev_path,
                   '-S', self.p.map_sam]
        # We have to tell bowtie2 if they we have FASTA files instead of FASTQ #
        if self.sample.format == 'fasta': options += ['-f']
        # Do the mapping #
        sh.bowtie2(*options)
        ## Create bam file, then sort it and finally index bamfile #
        sh.samtools('view', '-bt', self.contigs_fasta + '.fai', self.p.map_sam, '-o', self.p.map_bam)
        sh.samtools('sort', self.p.map_bam, self.p.map_s_bam.prefix_path)
        sh.samtools('index', self.p.map_s_bam)
        # Remove PCR duplicates #
        self.remove_duplicates()
        # Sort and index bam without duplicates #
        sh.samtools('sort', self.p.map_smd_bam, self.p.map_smds_bam.prefix_path)
        sh.samtools('index', self.p.map_smds_bam)
        # Compute coverage #
        sh.genomeCoverageBed('-ibam', self.p.map_smds_bam, _out=str(self.p.map_smds_coverage))
        # Clean up the ones we don't need #
        os.remove(self.p.map_sam)
        os.remove(self.p.map_bam)
        os.remove(self.p.map_smd_bam)

    def remove_duplicates(self):
        """Remove PCR duplicates with MarkDuplicates."""
        # Estimate size #
        mem_size = "16"
        # Run the command #
        sh.java('-Xmx%sg' % mem_size,
                '-XX:ParallelGCThreads=%s' % num_processors,
                '-XX:+CMSClassUnloadingEnabled',
                '-jar', gefes.repos_dir + 'bin/MarkDuplicates.jar',
                'INPUT=%s' % self.p.map_s_bam,
                'OUTPUT=%s' % self.p.map_smd_bam,
                'METRICS_FILE=%s' % self.p.map_smd_metrics,
                'AS=TRUE',
                'VALIDATION_STRINGENCY=LENIENT',
                'MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000',
                'REMOVE_DUPLICATES=TRUE')

    @property_cached
    def results(self):
        results = BowtieResults(self)
        if not results: raise Exception("You can't access results from Bowtie before running the mapping.")
        return results

###############################################################################
class BowtieResults(object):

    def __nonzero__(self): return self.p.coverage.exists
    def __init__(self, bowtie):
        self.bowtie = bowtie
        self.assembly = bowtie.assembly
        self.p = bowtie.p

    @property_pickled
    def statistics(self):
        """Uses the BEDTools genomeCoverageBed histogram output to determine mean
        coverage and percentage covered for each contig. Returns a dict with contig names as
        keys and containing 'percent_covered' with 'cov_mean' information inside.
        http://bedtools.readthedocs.org/en/latest/content/tools/genomecov.html"""
        headers = ['name', 'depth', 'count', 'length', 'fraction']
        frame = pandas.io.parsers.read_csv(self.p.map_smds_coverage, sep='\t', index_col=0, names=headers)
        out_dict = {}
        for contig in self.assembly.results.contigs:
            if contig.name not in frame.index:
                cov_mean = 0.0
                percent_covered = 0.0
            else:
                subframe = frame.loc[contig.name]
                assert round(subframe['fraction'].sum(), 4) == 1.0
                cov_mean = (subframe['depth'] * subframe['fraction']).sum()
                not_covered = subframe['depth'] == 0
                if not_covered.any(): percent_covered = float(100 - subframe.loc[not_covered]['fraction'])
                else:                 percent_covered = 100.0
            out_dict[contig.name] = {"cov_mean": cov_mean, "percent_covered": percent_covered}
        return out_dict

    @property_cached
    def coverage_mean(self):
        """The percentage covered in every contig. Dict with contig names as keys"""
        return {contig.name: self.statistics[contig.name]['cov_mean'] for contig in self.assembly.results.contigs}

    @property_cached
    def covered_fraction(self):
        """The mean coverage in every contig. Dict with contig names as keys"""
        return {contig.name: self.statistics[contig.name]['percent_covered'] for contig in self.assembly.results.contigs}