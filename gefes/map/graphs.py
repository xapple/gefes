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
    sep    = 'y'
    left   = 0.1
    title  = 'Distribution of mean coverages'
    x_label = 'Mean coverage of a contig'
    y_label = 'Number of contigs with this much mean coverage'

    def plot(self, **kwargs):
        # Data #
        counts = self.parent.coverage_mean.values()
        # Plot #
        fig = pyplot.figure()
        pyplot.hist(counts, bins=100, color='gray')
        axes = pyplot.gca()
        # Save it #
        self.save_plot(fig, axes, **kwargs)
        pyplot.close(fig)

################################################################################
class PercentCovered(Graph):
    """percent_covered"""
    short_name = 'percent_covered'
    sep    = 'y'
    left   = 0.1
    title  = 'Distribution of percentages covered'
    x_label = 'Percentage covered of a contig'
    y_label = 'Number of contigs with this much percent covered'

    def plot(self, **kwargs):
        # Data #
        counts = self.parent.covered_fraction.values()
        # Plot #
        fig = pyplot.figure()
        pyplot.hist(counts, bins=100, color='gray')
        axes = pyplot.gca()
        # Save it #
        self.save_plot(fig, axes, **kwargs)
        pyplot.close(fig)