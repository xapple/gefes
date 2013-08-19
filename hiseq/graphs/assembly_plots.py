# Built-in modules #

# Internal modules #
from hiseq.graphs import Graph

# Third party modules #
import pandas
from matplotlib import pyplot
from Bio import SeqIO

# Constants #
__all__ = ['ContigDist']

################################################################################
class ContigDist(Graph):
    """General distribution of the contigs for a pool"""
    short_name = 'contig_dist'

    def plot(self):
        # Data #
        lengths = map(len, SeqIO.parse(self.parent.p.amos_dir + 'contigs.fasta' , 'fasta'))
        values = pandas.Series(lengths)
        # Plot #
        fig = pyplot.figure()
        axes = values.hist(color='gray', bins=max(values))
        fig = pyplot.gcf()
        title = 'Distribution of contigs lengths after scafoloding for sample "%s"' % self.parent.pool.long_name
        axes.set_title(title)
        axes.set_xlabel('Number of nucleotides in sequence')
        axes.set_ylabel('Number of sequences with this length')
        axes.xaxis.grid(False)
        # Save it #
        self.save_plot(fig, axes, sep=('x'))