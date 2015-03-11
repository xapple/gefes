# Built-in modules #

# Internal modules #
from plumbing.graphs import Graph

# Third party modules #
from matplotlib import pyplot

# Constants #
__all__ = ['BinContigDistribution', 'BinNucleotideDistribution']

################################################################################
class BinContigDistribution(Graph):
    """bins_contig_dist"""
    short_name = 'bins_contig_dist'

    def plot(self, x_log=False, y_log=False, bins=250):
        # Data #
        counts = map(len, self.parent.bin_id_to_contig_ids.values())
        # Plot #
        fig = pyplot.figure()
        pyplot.hist(counts, bins=bins, color='gray')
        axes = pyplot.gca()
        # Information #
        title = 'Distribution of number of contigs in the bins'
        axes.set_title(title)
        axes.set_xlabel('Number of contigs in a bin')
        axes.set_ylabel('Number of bins with that many contigs in them')
        axes.xaxis.grid(False)
        # Add logarithm to axes #
        if x_log: axes.set_xscale('symlog')
        if y_log: axes.set_yscale('symlog')
        # Save it #
        self.save_plot(fig, axes, sep=('y'))
        # For convenience #
        return self

################################################################################
class BinNucleotideDistribution(Graph):
    """bins_nucleotide_dist"""
    short_name = 'bins_nucleotide_dist'

    def plot(self, x_log=False, y_log=True, bins=250):
        # Data #
        counts = [sum(map(len, b.contigs)) for b in self.parent.bins]
        # Plot #
        fig = pyplot.figure()
        pyplot.hist(counts, bins=bins, color='gray')
        axes = pyplot.gca()
        # Information #
        title = 'Distribution of the total nucleotide count in the bins'
        axes.set_title(title)
        axes.set_xlabel('Number of nucleotides in a bin')
        axes.set_ylabel('Number of bins with that many nucleotides in them')
        axes.xaxis.grid(False)
        # Add logarithm to axes #
        if x_log: axes.set_xscale('symlog')
        if y_log: axes.set_yscale('symlog')
        # Save it #
        self.save_plot(fig, axes, sep=('y'))
        # For convenience #
        return self