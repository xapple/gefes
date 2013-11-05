# Futures #
from __future__ import division

# Built-in modules #

# Internal modules #
from gefes.common.autopaths import AutoPaths
from gefes.fasta.paired import PairedFASTQ
from gefes.fasta.single import FASTQ

# Third party modules #
import sh

###############################################################################
class Cleaner(object):
    """Takes care of cleaning the raw reads."""

    all_paths = """
    /cleaned_fwd.fastq
    /cleaned_rev.fastq
    /cleaned_single.fastq
    /report.txt
    """

    def __repr__(self): return '<%s object of %s>' % (self.__class__.__name__, self.parent)
    def __len__(self): return len(self.pair)

    def __init__(self, pool):
        # Save parent #
        self.parent, self.pool = pool, pool
        # Auto paths #
        self.base_dir = self.parent.p.clean_dir
        self.p = AutoPaths(self.base_dir, self.all_paths)
        # Convenience objects #
        self.fwd = FASTQ(self.p.fwd)
        self.rev = FASTQ(self.p.rev)
        self.single = FASTQ(self.p.single)
        self.pair = PairedFASTQ(self.p.fwd, self.p.rev)

    def clean(self):
        # Cleanup #
        self.fwd.remove()
        self.rev.remove()
        self.single.remove()
        self.p.report.remove()
        # Call sickle #
        stats = sh.sickle("pe", "-f", self.pool.fwd, "-r", self.pool.rev, "-t", "sanger", "-o",
                          self.fwd, "-p", self.rev, "-s", self.single)
        # Write the report #
        with open(self.p.report, 'w') as handle: handle.write(str(stats))
        # Make a sanity check #
        assert self.kept + self.discarded == len(self.parent)

    @property
    def stats(self):
        # Parse the report file #
        return self.p.report.contents

    @property
    def kept(self):
        return len(self)

    @property
    def discarded(self):
        return len(self.parent) - len(self)