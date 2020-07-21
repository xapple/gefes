# Futures #
from __future__ import division

# Built-in modules #
import itertools

# Internal modules #
from gefes.preprocess import QualityChecker, QualityResults

# First party modules #
from plumbing.common import moving_average
from plumbing.cache import property_cached

# Third party modules #

###############################################################################
class SlidingWindow(QualityChecker):
    """
    Takes care of checking the PHRED score of the raw reads
    and will discard or trim bad ones.

    Other interesting approach here:

    https://github.com/novigit/broCode/blob/master/pbamp/PacBioTrimmer.py#L29
    """

    window_size = 10 # Size of the window we will slide along the read
    threshold   = 20 # This is a PHRED score threshold
    min_length  = 71 # Minimum number of remaining base pairs
    discard_N   = True

    def run(self):
        # Count #
        self.discarded = 0
        # Do it #
        with self.dest as output, self.singletons as singles:
            for read_pair in self.source.progress:
                one = self.trim_read(read_pair[0])
                two = self.trim_read(read_pair[1])
                if one and two: output.add_pair((one, two))
                elif one: singles.add_seq(one)
                elif two: singles.add_seq(two)
                else: self.discarded += 1
        # Make sanity checks #
        assert len(self.source) == len(self.dest) + len(self.singletons) + self.discarded
        # Return result #
        return self.results

    def trim_read(self, read):
        # First we remove base pairs strictly below the threshold on both sides #
        phred = read.letter_annotations["phred_quality"]
        above_yes_no = [True if x > self.threshold else False for x in phred]
        if True not in above_yes_no: return None
        new_start = above_yes_no.index(True)
        new_end = list(reversed(above_yes_no)).index(True)
        if new_end == 0: read = read[new_start:]
        if new_end != 0: read = read[new_start:-new_end]
        if len(read.seq) < self.min_length: return None
        # Now we run our moving average #
        phred = read.letter_annotations["phred_quality"]
        averaged = moving_average(phred, self.window_size, 'copy_padding_and_cut')
        # And then search for the longest stretch above the threshold #
        stretches = itertools.groupby(averaged, lambda x: x>self.threshold)
        stretches = [Stretch(above, scores) for above, scores in stretches]
        if not stretches: return None
        # This is a bit clumsy but we need to calculate the starts and ends #
        start = 0
        for s in stretches:
            s.start = start
            start = start + len(s)
            s.end = start
            s.seq = read[s.start:s.end]
        # Remove the below ones #
        stretches = [s for s in stretches if s.above]
        if not stretches: return None
        # Discard N letters #
        if self.discard_N: stretches = [s for s in stretches if 'N' not in s.seq]
        if not stretches: return None
        # Get the longest one #
        longest = max(stretches, key = lambda x: len(x))
        # Check the length #
        if len(longest.seq) < self.min_length: return None
        else: return longest.seq

    @property_cached
    def results(self):
        results = SlidingWindowResults(self, self.source, self.dest, self.singletons)
        if not results:
            raise Exception("You can't access results from the quality check before running the algorithm.")
        return results

###############################################################################
class SlidingWindowResults(QualityResults): pass

###############################################################################
class Stretch(object):
    """
    An interval within a DNA sequence. Here, a continuous stretch where
    every base is above a given threshold.
    """

    def __repr__(self): return "%s from %i to %i" % (self.above, self.start, self.end)
    def __len__(self): return len(self.scores)
    def __init__(self, above, scores):
        self.above  = above
        self.scores = list(scores)
        self.start  = -1
        self.end    = -1

###############################################################################
def test():
    """Test that our implementation is correct on one sequence"""
    # Object #
    checker = QualityChecker(None, None)
    checker.window_size = 10
    checker.threshold   = 21
    checker.min_length  = 5
    checker.discard_N   = True
    # Dummy test sequence #
    scores = "10 11 13 22 23 24 10 10 9 8 7 9 8 9 5 2 5 8 9 8 9 30 33 30 31 32 33 31 33 33 31 33 32 33 32 32 33 32 2 3 2 3 2 1 3 2 1 23 23 23 10 10 9 9"
    seq    = "A  T  C  G  T  T  G  A  C G G A G T G T A A C T C G  A  T  G  A  C  T  T  G  T  C  A  A  C  T  G  G  T A G G G T C A A C  T  G  A  T  C A"
    scores = map(int, scores.split())
    seq    = ''.join(seq.split())
    assert len(seq) == len(scores)
    # Make them into biopython objects #
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
