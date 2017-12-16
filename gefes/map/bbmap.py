# Futures #
from __future__ import division

# Built-in modules #
import os, sys

# Internal modules #
from gefes.map import Mapper, MapperResults

# First party modules #
from plumbing.cache import property_cached
from plumbing.slurm import num_processors

# Third party modules #
import sh

# Constants #
home = os.environ.get('HOME', '~') + '/'

###############################################################################
class BBMap(Mapper):
    """Uses https://sourceforge.net/projects/bbmap/
    Example: bbmap.sh in1=reads1.fq in2=reads2.fq out=mapped.sam ref=ecoli.fa nodisk"""

    short_name = 'bbmap'
    long_name  = 'BBMap_37.76'
    executable = 'bbmap.sh'
    url        = 'https://sourceforge.net/projects/bbmap/'
    dependencies = ['java']

    def pre_run(self, verbose=True):
        """Check an index exist."""
        # FAI Index #
        if not os.path.exists(self.contigs_fasta + '.fai'):
            if verbose: print "Making samtools index"; sys.stdout.flush()
            self.contigs_fasta.index_samtools()

    def run(self, verbose=True, cpus=None):
        # Variable threads #
        if cpus is None: cpus = num_processors
        # Make our options #
        options = ['ref=%s' % self.assembly.results.contigs_fasta,
                   'in1=%s' % self.sample.clean.fwd,
                   'in2=%s' % self.sample.clean.rev,
                   'out=%s' % self.p.map_sam,
                   'threads=%i' % cpus,
                   'nodisk']
        # Messages #
        if verbose: print "Launching BBMap on sample '%s'..." % self.sample.name
        if verbose: print "Mapping against assembly '%s'." % self.assembly
        sys.stdout.flush()
        # Do the mapping #
        bbmap = sh.Command(home + "programs/bbmap/bbmap.sh")
        bbmap(*options)

    @property_cached
    def results(self):
        results = BBMapResults(self)
        if not results: raise Exception("You can't access results from BBMap before running the mapping.")
        return results

###############################################################################
class BBMapResults(MapperResults): pass