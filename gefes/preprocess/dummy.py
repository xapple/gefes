# Futures #
from __future__ import division

# Built-in modules #

# Internal modules #
from gefes.preprocess import QualityChecker, QualityResults

# First party modules #
from plumbing.cache import property_cached

# Third party modules #

###############################################################################
class DummyCleaner(QualityChecker):
    """For when it's already been done by someone else."""

    @property_cached
    def results(self):
        results = DummyCleaner(self, self.source, self.dest, self.singletons)
        return results

###############################################################################
class DummyCleanerResults(QualityResults):
    pass