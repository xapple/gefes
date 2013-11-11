# Built-in modules #

# Internal modules #
from gefes.graphs import Graph

# Third party modules #
import pandas
from matplotlib import pyplot

# Constants #
__all__ = ['ContigDist']

################################################################################
class ContigDist(Graph):
    """General distribution of the contigs for a pool"""
    short_name = 'contig_dist'

    def plot(self):
        # Data #
        values = pandas.Series(self.parent.contigs_fasta.lengths)
        # Plot #
        fig = pyplot.figure()
        axes = values.hist(color='gray', bins=2000)
        fig = pyplot.gcf()
        title = 'Distribution of contigs lengths after scaffolding for sample "%s"' % self.parent.pool.long_name
        axes.set_title(title)
        axes.set_xlabel('Number of nucleotides in sequence')
        axes.set_ylabel('Number of sequences with this length')
        axes.xaxis.grid(False)
        # Save it #
        self.save_plot(fig, axes, sep=('x'))