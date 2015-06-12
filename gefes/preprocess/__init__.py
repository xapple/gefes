# Futures #
from __future__ import division

# Built-in modules #

# Internal modules #
from fasta import PairedFASTQ, FASTQ

# Third party modules #

###############################################################################
class QualityChecker(object):
    """The base QualityChecker others should inherit from."""

    def __repr__(self): return '<%s object of %s>' % (self.__class__.__name__, self.source)
    def __len__(self): return len(self.pair)

    def __init__(self, source, dest=None):
        # Basic #
        self.source = source
        self.dest = dest
        # Default case #
        if dest is None:
            self.dest = PairedFASTQ(self.source.fwd.prefix_path + '_clean.fastq',
                                    self.source.rev.prefix_path + '_clean.fastq')
        # Single read #
        self.singletons = FASTQ(self.dest.fwd.directory + 'singletons.fastq')

    def run(self): raise NotImplementedError('You have to implement this method.')

    @property
    def results(self): raise NotImplementedError('You have to implement this method.')

###############################################################################
class QualityResults(object):
    """The base QualityResults others should inherit from."""

    all_paths = """
    /lorem
    """

    def __nonzero__(self): return bool(self.dest)

    def __init__(self, checker, source, dest, singletons):
        self.checker    = checker
        self.source     = source
        self.dest       = dest
        self.singletons = singletons

    @property
    def ratio_kept(self):
        return len(self.dest) / len(self.source)