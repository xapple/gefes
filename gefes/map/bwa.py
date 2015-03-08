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
from shell_command import shell_output

###############################################################################
class Bwa(Mapper):
    """Use BWA at http://bio-bwa.sourceforge.net
    to maps reads from a Sample object back to the contigs of an
    Assembly object. Expects version 0.7.12-r1039 of bwa.
    """

    short_name = 'bwa-mem'
    long_name  = 'BWA-MEM 0.7.12-r1039'
    executable = 'bwa'

    def run(self, verbose=True):
        # Make our options #
        options = ['mem',
                   self.assembly.results.contigs_fasta,
                   self.sample.clean.fwd,
                   self.sample.clean.rev,
                   '-t', num_processors]
        # Do the mapping #
        if verbose: print "Launching BWA on sample '%s'..." % self.sample.name ; sys.stdout.flush()
        shell_output('bwa' + ' '.join(options) + " > " + self.p.map_sam)

    @property_cached
    def results(self):
        results = BwaResults(self)
        if not results: raise Exception("You can't access results from Bowtie before running the mapping.")
        return results

###############################################################################
class BwaResults(MapperResults): pass