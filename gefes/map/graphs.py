# Built-in modules #
from collections import Counter

# Internal modules #
from plumbing.graphs import Graph

# Third party modules #
from matplotlib import pyplot

# Constants #
__all__ = ['MeanCoverage', 'PercentCovered']

################################################################################
class MeanCoverage(Graph):
    """mean_coverage"""
    short_name = 'mean_coverage'
    sep = ('x')

    def plot(self, **kwargs):
        # Data #
        counts = self.parent.coverage_mean.values()
        # Plot #
        fig = pyplot.figure()
        pyplot.hist(counts, bins=100, color='gray')
        axes = pyplot.gca()
        # Information #
        title = 'Distribution of mean coverages'
        axes.set_title(title)
        axes.set_xlabel('Mean coverage of a contig')
        axes.set_ylabel('Number of contigs with this much mean coverage')
        # Save it #
        self.save_plot(fig, axes, **kwargs)
        pyplot.close(fig)

class PercentCovered(Graph):
    """percent_covered"""
    short_name = 'percent_covered'
    sep = ('x')

    def plot(self, **kwargs):
        # Data #
        counts = self.parent.covered_fraction.values()
        # Plot #
        fig = pyplot.figure()
        pyplot.hist(counts, bins=100, color='gray')
        axes = pyplot.gca()
        # Information #
        title = 'Distribution of percentages covered'
        axes.set_title(title)
        axes.set_xlabel('Percentage covered of a contig')
        axes.set_ylabel('Number of contigs with this much percent covered')
        # Save it #
        self.save_plot(fig, axes, **kwargs)
        pyplot.close(fig)