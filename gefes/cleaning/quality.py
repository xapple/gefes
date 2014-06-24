# Built-in modules #

# Internal modules #
from gefes.common import moving_average
from gefes.common.autopaths import AutoPaths
from gefes.fasta.paired import PairedFASTQ
from gefes.fasta.single import FASTQ

# Third party modules #

###############################################################################
class QualityChecker(object):
    """Takes care of checking the PHRED score of the raw reads and discarding bad ones."""

    window_size = 10
    qual_threshold = 20

    all_paths = """
    /cleaned_fwd.fastq
    /cleaned_rev.fastq
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
        self.paired = PairedFASTQ(self.p.fwd, self.p.rev)

    def run(self, discard_N=False):
        # Cleanup #
        self.fwd.remove()
        self.rev.remove()
        # Call sickle #
        self.paired.create()
        for read_pair in self.pool.pair:
            if discard_N:
                if "N" in read_pair[0] or "N" in read_pair[1]: continue
            one = moving_average(read_pair[0].letter_annotations["phred_quality"], self.window_size)
            two = moving_average(read_pair[1].letter_annotations["phred_quality"], self.window_size)
            if any([value < self.qual_threshold for value in one]): continue
            if any([value < self.qual_threshold for value in two]): continue
            self.paired.add_pair(read_pair)
        self.paired.close()
        # Make sanity checks #
        assert len(self.fwd) == len(self.rev)