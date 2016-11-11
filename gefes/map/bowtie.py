# Futures #
from __future__ import division

# Built-in modules #
import sys

# Internal modules #
from gefes.map import Mapper, MapperResults

# First party modules #
from plumbing.cache import property_cached
from plumbing.slurm import num_processors

# Third party modules #
import sh

###############################################################################
class Bowtie(Mapper):
    """Uses Bowtie2 to maps reads from a Sample object back to an Assembly object.
    Expects version 2.2.5.
    By default, if a read can map in several spots, from the manual:
    'The best alignment found is reported (randomly selected from among best if tied)'
    """

    short_name = 'bowtie'
    long_name  = 'Bowtie2 2.2.5'
    executable = 'bowtie2'
    url        = 'http://bowtie-bio.sourceforge.net/bowtie2/index.shtml'
    dependencies = []

    def run(self, verbose=True, cpus=None):
        # Check both type of indexes exist #
        self.pre_run()
        # Variable threads #
        if cpus is None: cpus = num_processors
        # Make our options #
        self.options = ['-p', str(cpus),
                        '-x', self.assembly.results.contigs_fasta,
                        '-1', self.sample.clean.fwd,
                        '-2', self.sample.clean.rev,
                        '-S', self.p.map_sam]
        # We have to tell Bowtie2 if they we have FASTA files instead of FASTQ #
        if self.sample.clean.format == 'fasta': self.options += ['-f']
        # Do the mapping #
        if verbose: print "Launching Bowtie on sample '%s' with %i cores" % (self.sample.name, cpus)
        if verbose: print "Mapping against assembly '%s'." % self.assembly
        sys.stdout.flush()
        sh.bowtie2(*self.options, _out=self.p.stdout.path, _err=self.p.stderr.path)
        # Create bam file, then sort it and finally index the bamfile #
        self.post_run(cpus=cpus)

    @property_cached
    def results(self):
        results = BowtieResults(self)
        if not results: raise Exception("You can't access results from Bowtie before running the mapping.")
        return results

###############################################################################
class BowtieResults(MapperResults): pass