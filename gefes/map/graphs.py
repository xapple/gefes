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

    def plot(self, x_log=False, y_log=False):
        # Data #
        counts = self.parent.coverage_mean.values()
        # Plot #
        fig = pyplot.figure()
        pyplot.hist(counts, bins=100, color='k')
        axes = pyplot.gca()
        # Information #
        title = 'Distribution of mean coverages'
        axes.set_title(title)
        axes.set_xlabel('Mean coverage of a contig')
        axes.set_ylabel('Number of contigs with this much mean coverage')
        axes.xaxis.grid(False)
        # Add logarithm to axes #
        if x_log: axes.set_xscale('symlog')
        if y_log: axes.set_yscale('symlog')
        # Save it #
        self.save_plot(fig, axes, sep=('x'))

class PercentCovered(Graph):
    """percent_covered"""
    short_name = 'percent_covered'

    def plot(self, x_log=False, y_log=False):
        # Data #
        counts = self.parent.covered_fraction.values()
        # Plot #
        fig = pyplot.figure()
        pyplot.hist(counts, bins=100, color='k')
        axes = pyplot.gca()
        # Information #
        title = 'Distribution of percentages covered'
        axes.set_title(title)
        axes.set_xlabel('Percentage covered of a contig')
        axes.set_ylabel('Number of contigs with this much percent covered')
        axes.xaxis.grid(False)
        # Add logarithm to axes #
        if x_log: axes.set_xscale('symlog')
        if y_log: axes.set_yscale('symlog')
        # Save it #
        self.save_plot(fig, axes, sep=('x'))