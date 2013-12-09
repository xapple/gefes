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
    """General distribution of the contig lengths"""
    short_name = 'contig_length'

    def plot(self):
        # Data #
        values = pandas.Series(self.parent.contigs_fasta.lengths)
        # Plot #
        fig = pyplot.figure()
        axes = values.hist(color='gray', bins=max(values))
        fig = pyplot.gcf()
        title = 'Distribution of contig lengths for %i contigs' % self.parent.contigs_fasta.count
        axes.set_title(title)
        axes.set_xlabel('Number of nucleotides in sequence')
        axes.set_ylabel('Number of sequences with this length')
        axes.xaxis.grid(False)
        axes.set_yscale('symlog')
        axes.set_xscale('symlog')
        # Save it #
        self.save_plot(fig, axes, sep=('x'))