# Built-in modules #
import re

# Internal modules #
from gefes.preprocess import QualityChecker, QualityResults

# First party modules #
from plumbing.cache import property_cached

# Third party modules #
import sh

###############################################################################
class Sickle(QualityChecker):
    """Takes care of running the sickle program that removes low quality reads.
    The 'singletons' file contains reads that passed the filter in either
    the forward or reverse direction, but not the other"""

    def run(self):
        # Cleanup #
        self.fwd.remove()
        self.rev.remove()
        self.single.remove()
        self.p.report.remove()
        # Call sickle #
        sh.sickle("pe", "-f", self.pool.fwd, "-r", self.pool.rev, "-t", "sanger", "-n",
                  "-o", self.fwd, "-p", self.rev, "-s", self.single, _out=str(self.p.report))
        # Make sanity checks #
        self.parse_stats()
        assert len(self.fwd) == len(self.rev)
        assert self.paired_records_kept == len(self.fwd)
        assert self.single_records_kept == len(self.single)
        assert self.kept + self.discarded == len(self.pool.fwd)
        # Make sanity checks #
        assert len(self.source) == self.discarded + len(self.singletons) + len(self.dest)
        # Return result #
        return self.results

    @property_cached
    def results(self):
        results = SickleResults(self.source, self.dest, self.singletons)
        if not results: raise Exception("You can't access results from sickle before running the algorithm.")
        return results

###############################################################################
class SickleResults(QualityResults):

    @property_cached
    def stats(self):
        """Parse the report file for statistics"""
        result = {}
        result['paired_records_kept']      = int(re.findall('^FastQ paired records kept (.+) .+$',      self.p.report.contents, re.M))
        result['single_records_kept']      = int(re.findall('^FastQ single records kept (.+) .+$',      self.p.report.contents, re.M))
        result['paired_records_discarded'] = int(re.findall('^FastQ paired records discarded (.+) .+$', self.p.report.contents, re.M))
        result['single_records_discarded'] = int(re.findall('^FastQ single records discarded (.+) .+$', self.p.report.contents, re.M))
        return result
