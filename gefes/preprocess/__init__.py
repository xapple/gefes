# Futures #
from __future__ import division

# Built-in modules #

# Internal modules #

# First party #
from fasta import PairedFASTQ, FASTQ
from plumbing.autopaths import AutoPaths

# Third party modules #

###############################################################################
class QualityChecker(object):
    """The base QualityChecker others should inherit from."""

    all_paths = """
    /report.txt
    """

    def __repr__(self): return '<%s object on %s>' % (self.__class__.__name__, self.source)
    def __len__(self): return len(self.pair)

    def __init__(self, result_dir, source, dest=None, singletons=None):
        # Basic #
        self.result_dir = result_dir
        self.source     = source
        self.dest       = dest
        self.singletons = singletons
        # Destination #
        if self.dest is None:
            self.dest = PairedFASTQ(self.source.fwd.prefix_path + '_clean.fastq',
                                    self.source.rev.prefix_path + '_clean.fastq')
        # Single read #
        if self.singletons is None:
            self.singletons = FASTQ(self.dest.fwd.directory + 'singletons.fastq')
        # Auto paths #
        self.base_dir = self.result_dir # + self.short_name + '/'
        self.p = AutoPaths(self.base_dir, self.all_paths)

    def run(self): raise NotImplementedError('You have to implement this method.')

    @property
    def results(self): raise NotImplementedError('You have to implement this method.')

###############################################################################
class QualityResults(object):
    """The base QualityResults others should inherit from."""

    def __nonzero__(self): return bool(self.dest)

    def __init__(self, checker, source, dest, singletons):
        self.checker    = checker
        self.source     = source
        self.dest       = dest
        self.singletons = singletons

    @property
    def ratio_kept(self):
        return len(self.dest) / len(self.source)