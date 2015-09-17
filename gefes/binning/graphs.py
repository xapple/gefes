# Built-in modules #

# Internal modules #
from plumbing.graphs import Graph

# Third party modules #
from matplotlib import pyplot

# Constants #
__all__ = ['BinContigDistribution', 'BinNucleotideDistribution', 'BinGenesPredictedDist']

################################################################################
class BinContigDistribution(Graph):
    """bins_contig_dist"""
    short_name = 'bins_contig_dist'
    sep = ('y')

    def plot(self, bins=250, **kwargs):
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
        self.save_plot(fig, axes, **kwargs)
        pyplot.close(fig)
        # For convenience #
        return self

################################################################################
class BinNucleotideDistribution(Graph):
    """bins_nucleotide_dist"""
    short_name = 'bins_nucleotide_dist'
    sep = ('y')
    y_scale = 'symlog'

    def plot(self, bins=250, **kwargs):
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
        # Save it #
        self.save_plot(fig, axes, **kwargs)
        pyplot.close(fig)
        # For convenience #
        return self

################################################################################
class BinGenesPredictedDist(Graph):
    """bins_genes_predicted_dist"""
    short_name = 'bins_genes_predicted_dist'
    sep = ('y')
    y_scale = 'symlog'

    def plot(self, bins=250, **kwargs):
        # Data #
        counts = [b.evaluation.results[''] for b in self.parent.bins]
        # Plot #
        fig = pyplot.figure()
        pyplot.hist(counts, bins=bins, color='gray')
        axes = pyplot.gca()
        # Information #
        title = 'Distribution of the total nucleotide count in the bins'
        axes.set_title(title)
        axes.set_xlabel('Number of nucleotides in a bin')
        axes.set_ylabel('Number of bins with that many nucleotides in them')
        # Save it #
        self.save_plot(fig, axes, **kwargs)
        pyplot.close(fig)
        # For convenience #
        return self