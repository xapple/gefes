# Built-in modules #

# Internal modules #
from gefes.graphs import Graph

# Third party modules #
import pandas
from matplotlib import pyplot

# Constants #
__all__ = ['CleanLengthDist']

################################################################################
class CleanLengthDist(Graph):
    """The distribution of the reads after the cleaning"""
    short_name = 'clean_length_dist'

    def plot(self):
        # Data #
        values = pandas.Series(len(s) for s in self.parent.cleaner.fwd)
        # Plot #
        fig = pyplot.figure()
        axes = values.hist(color='gray', bins=100)
        fig = pyplot.gcf()
        title = 'Distribution of sequence lengths after cleaning for pool "%s"' % self.parent.long_name
        axes.set_title(title)
        axes.set_xlabel('Number of nucleotides in sequence')
        axes.set_ylabel('Number of sequences with this length')
        axes.xaxis.grid(False)
        # Save it #
        self.save_plot(fig, axes, sep=('x'))