# Futures #
from __future__ import division

# Built-in modules #

# Internal modules #
from gefes.map import Mapper, MapperResults

# First party modules #
from plumbing.cache import property_cached
from plumbing.slurm import num_processors

# Third party modules #
import sh

###############################################################################
class Bowtie(Mapper):
    """Uses Bowtie2 at http://bowtie-bio.sourceforge.net/bowtie2/index.shtml
    to maps reads from a Sample object back to an Assembly object.
    Expects version 2.2.4.
    """

    short_name = 'bowtie'
    long_name  = 'Bowtie2 2.2.4'
    executable = 'bowtie2'

    def run(self, verbose=True):
        # Check both type of indexes exist #
        self.pre_run()
        # Make our options #
        options = ['-p', num_processors,
                   '-x', self.assembly.results.contigs_fasta,
                   '-1', self.sample.clean.fwd,
                   '-2', self.sample.clean.rev,
                   '-S', self.p.map_sam]
        # We have to tell bowtie2 if they we have FASTA files instead of FASTQ #
        if self.sample.format == 'fasta': options += ['-f']
        # Do the mapping #
        if verbose: print "Launching Bowtie on sample '%s'..." % self.sample.name
        sh.bowtie2(*options)
        # Create bam file, then sort it and finally index the bamfile #
        self.post_run()

    @property_cached
    def results(self):
        results = BowtieResults(self)
        if not results: raise Exception("You can't access results from Bowtie before running the mapping.")
        return results

###############################################################################
class BowtieResults(MapperResults): pass