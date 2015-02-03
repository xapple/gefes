# Built-in modules #

# Internal modules #
from plumbing.graphs import Graph

# Third party modules #
from matplotlib import pyplot

# Constants #
__all__ = ['BinSizeDistribution']

################################################################################
class BinSizeDistribution(Graph):
    """bins_size_dist"""
    short_name = 'bins_size_dist'

    def plot(self, x_log=False, y_log=True):
        # Data #
        counts = map(len,self.parent.bin_id_to_contig_ids.values())
        # Plot #
        fig = pyplot.figure()
        pyplot.hist(counts, bins=100, color='gray')
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