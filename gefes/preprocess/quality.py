# Built-in modules #
import itertools

# Internal modules #
from plumbing.common import moving_average

# Third party modules #

###############################################################################
class QualityChecker(object):
    """Takes care of checking the PHRED score of the raw reads
    and will discard or trim bad ones."""

    window_size = 10
    threshold = 20
    min_length = 50
    discard_N = True

    def __repr__(self): return '<%s object of %s>' % (self.__class__.__name__, self.parent)
    def __len__(self): return len(self.pair)

    def __init__(self, source, dest):
        self.source = source
        self.dest = dest

    def run(self):
        # Do it #
        with self.dest as output:
            for read_pair in self.source:
                one = self.trim_read(read_pair[0])
                two = self.trim_read(read_pair[1])
                output.add_pair((one, two))
        # Make sanity checks #
        assert len(self.source) == len(self.dest)

    def trim_read(self, read):
        # First we remove base pairs strictly below the threshold on both sides #
        phred = read.letter_annotations["phred_quality"]
        above_yes_no = [True if x > self.threshold else False for x in phred]
        new_start = above_yes_no.index(True)
        new_end = list(reversed(above_yes_no)).index(True)
        read = read[new_start:-new_end]
        phred = read.letter_annotations["phred_quality"]
        # Now we run our moving average #
        averaged = moving_average(phred, self.window_size, 'copy_padding_and_cut')
        # And then search for the longest stretch above the threshold #
        stretches = itertools.groupby(averaged, lambda x: x>self.threshold)
        stretches = [Stretch(above, scores) for above, scores in stretches]
        # This is a bit clumsy but we need to calculate the starts and ends #
        start = 0
        for s in stretches:
            s.start = start
            start = start + len(s)
            s.end = start
            s.seq = read[s.start:s.end]
        # Discard N letters #
        if self.discard_N: stretches = [s for s in stretches if 'N' not in s.seq]
        # Get the longest one #
        longest = max([s for s in stretches if s.above], key = lambda x: len(x))
        # Check the length #
        if len(longest.seq) < self.min_length: return None
        else: return longest.seq

###############################################################################
class Stretch(object):
    """An interval within a DNA sequence. Here, a continuous stretch where
    every base is above a given threshold"""

    def __repr__(self): return "%s from %i to %i" % (self.above, self.start, self.end)
    def __len__(self): return len(self.scores)
    def __init__(self, above, scores):
        self.above = above
        self.scores = list(scores)
        self.start = -1
        self.end = -1

###############################################################################
class QualityResults(object):

    all_paths = """
    /lorem
    """

    def __init__(self, base_dir):
        self.base_dir = base_dir

    @property
    def ratio_discarded(self):
        return 1 - (len(self.dest) / len(self.source))

###############################################################################
def test():
    # Object #
    checker = QualityChecker(None, None)
    checker.window_size = 10
    checker.threshold = 21
    checker.min_length = 5
    checker.discard_N = True
    # Dummy sequence #
    scores = "10 11 13 22 23 24 10 10 9 8 7 9 8 9 5 2 5 8 9 8 9 30 33 30 31 32 33 31 33 33 31 33 32 33 32 32 33 32 2 3 2 3 2 1 3 2 1 23 23 23 10 10 9 9"
    seq    = "A  T  C  G  T  T  G  A  C G G A G T G T A A C T C G  A  T  G  A  C  T  T  G  T  C  A  A  C  T  G  G  T A G G G T C A A C  T  G  A  T  C A"
    scores = map(int, scores.split())
    seq    = ''.join(seq.split())
    assert len(seq) == len(scores)
    # Make into biopython object #
    from Bio.SeqRecord import SeqRecord
    from Bio.Seq import Seq
    read = SeqRecord(Seq(seq), id="test", name="test", description="test")
    read.letter_annotations["phred_quality"] = scores
    # Trim it #
    trimmed = checker.trim_read(read)
    # Check result #
    correct = "30 33 30 31 32 33 31 33 33 31 33 32 33 32 32"
    correct = "G  A  T  G  A  C  T  T  G  T  C  A  A  C  T "
    correct = ''.join(correct.split())
    assert str(trimmed.seq) == correct
