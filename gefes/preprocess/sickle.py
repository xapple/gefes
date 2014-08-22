# Built-in modules #
import re

# Internal modules #

# Third party modules #
import sh

###############################################################################
class Sickle(object):
    """Takes care of running the sickle program that removes low quality reads.
    The "singles" file contains reads that passed filter in either the forward or reverse direction, but not the other"""

    def __repr__(self): return "<%s object on '%s'>" % (self.__class__.__name__, self.source)

    def __init__(self, source, dest):
        self.source = source
        self.dest = dest

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

    def parse_stats(self):
        # Parse the report file #
        self.paired_records_kept = int(re.findall('^FastQ paired records kept (.+) .+$', self.p.report.contents, re.M))
        self.single_records_kept = int(re.findall('^FastQ single records kept (.+) .+$', self.p.report.contents, re.M))
        self.paired_records_discarded = int(re.findall('^FastQ paired records discarded (.+) .+$', self.p.report.contents, re.M))
        self.single_records_discarded = int(re.findall('^FastQ single records discarded (.+) .+$', self.p.report.contents, re.M))

###############################################################################
class SickleResults(object):

    all_paths = """
    /cleaned_fwd.fastq
    /cleaned_rev.fastq
    /cleaned_single.fastq
    /report.txt
    """

    def __init__(self, base_dir):
        self.base_dir = base_dir

    @property
    def kept(self):
        return (self.paired_records_kept*2) + self.single_records_kept

    @property
    def discarded(self):
        return (self.paired_records_discarded*2) + self.single_records_discarded