# Futures #
from __future__ import division

# Built-in modules #

# Internal modules #
from gefes.graphs import Graph

# Third party modules #
import pandas, matplotlib
from matplotlib import pyplot

# Constants #
__all__ = ['ContigDist', 'ContigCumsum']

################################################################################
class ContigDist(Graph):
    """General distribution of the contig lengths"""
    short_name = 'contig_length'

    def plot(self):
        # Data #
        values = pandas.Series(self.parent.contigs_fasta.lengths)
        # Plot #
        fig = pyplot.figure()
        axes = values.hist(color='gray', bins=1000)
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

################################################################################
class ContigCumsum(Graph):
    """Cumulative sum of nuclotides by contig lengths (c.f. brouillon)"""
    short_name = 'contig_cumsum'

    def plot(self):
        # Data #
        values = pandas.Series(self.parent.contigs_fasta.lengths)
        values.sort()
        cumsum = values.cumsum()
        cumsum = cumsum / max(cumsum)
        # Plot #
        fig = pyplot.figure()
        pyplot.step(values, cumsum, '-k')
        title = 'Cumulative sum of contig lengths sorted from smallest to largest.'
        axes = pyplot.gca()
        axes.set_title(title)
        axes.set_xlabel('Length of sequence in nucleotides')
        axes.set_ylabel('Cumulative sum of nucleotides up to this length of contig')
        axes.set_xscale('symlog')
        axes.yaxis.grid(True)
        axes.xaxis.grid(True)
        # Percentage #
        def percentage(x, pos): return '%1.0f%%' % (x*100.0)
        axes.yaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(percentage))
        # Save it #
        self.save_plot(fig, axes, left=0.09)
        pyplot.close(fig)