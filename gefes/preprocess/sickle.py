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
    the forward or reverse direction, but not the other.
    Expects version 1.33
    """

    short_name = 'sickle'
    long_name  = 'Sickle cleaner v1.33'
    executable = 'sickle133'
    url        = 'https://github.com/najoshi/sickle/'
    dependencies = []

    def run(self):
        # Check version #
        assert "version 1.33" in sh.sickle133('--version')
        # Call sickle #
        sh.sickle133("pe", "-n",
                     "-f", self.source.fwd, "-r", self.source.rev,
                     "-o", self.dest.fwd,   "-p", self.dest.rev,
                     "-s", self.singletons,
                     "-t", "sanger",
                     _out=self.p.report.path)
        # Make sanity checks #
        assert len(self.dest.fwd) == len(self.dest.rev)
        assert len(self.source) == len(self.dest) + len(self.singletons) + self.discarded
        # Return result #
        return self.results

    @property_cached
    def results(self):
        results = SickleResults(self, self.source, self.dest, self.singletons)
        if not results: raise Exception("You can't access results from sickle before running the algorithm.")
        return results

###############################################################################
class SickleResults(QualityResults):

    @property_cached
    def stats(self):
        """Parse the report file for statistics"""
        patterns = {'paired_records_kept':      '^FastQ paired records kept (.+) .+$',
                    'single_records_kept':      '^FastQ single records kept (.+) .+$',
                    'paired_records_discarded': '^FastQ paired records discarded (.+) .+$',
                    'single_records_discarded': '^FastQ single records discarded (.+) .+$'}
        return {k: int(re.findall(v, self.p.report.contents, re.M)) for k,v in patterns.items()}