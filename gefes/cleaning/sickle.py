# Built-in modules #
import re

# Internal modules #
from plumbing.autopaths import AutoPaths
from fasta import PairedFASTQ
from fasta import FASTQ

# Third party modules #
import sh

###############################################################################
class Sickle(object):
    """Takes care of running the sickle program that removes low quality reads.
    The "singles" file contains reads that passed filter in either the forward or reverse direction, but not the other"""

    all_paths = """
    /cleaned_fwd.fastq
    /cleaned_rev.fastq
    /cleaned_single.fastq
    /report.txt
    """

    def __repr__(self): return '<%s object of %s>' % (self.__class__.__name__, self.parent)
    def __len__(self): return len(self.pair)

    def __init__(self, cleaner):
        # Save parent #
        self.parent, self.cleaner = cleaner, cleaner
        self.pool = self.cleaner.pool
        # Auto paths #
        self.base_dir = self.parent.base_dir
        self.p = AutoPaths(self.base_dir, self.all_paths)
        # Files #
        self.fwd = FASTQ(self.p.fwd)
        self.rev = FASTQ(self.p.rev)
        self.single = FASTQ(self.p.single)
        self.pair = PairedFASTQ(self.p.fwd, self.p.rev)

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

    @property
    def kept(self):
        return (self.paired_records_kept*2) + self.single_records_kept

    @property
    def discarded(self):
        return (self.paired_records_discarded*2) + self.single_records_discarded