# Built-in modules #

# Internal modules #
from plumbing.graphs import Graph

# Third party modules #
import numpy
from matplotlib import pyplot

# Constants #
__all__ = ['BinContigDistribution', 'BinNucleotideDistribution', 'BinGenesPredictedDist']

################################################################################
class BinContigDistribution(Graph):
    """Bin number of contigs distribution"""
    short_name = 'bins_contig_dist'
    sep        = 'x'
    y_scale    = 'symlog'

    def plot(self, bins=80, **kwargs):
        # Data #
        counts = map(len, self.parent.bin_id_to_contig_ids.values())
        # Linear bins in logarithmic space #
        if 'log' in kwargs.get('x_scale', ''):
            start, stop = numpy.log10(1), numpy.log10(max(counts))
            bins = list(numpy.logspace(start=start, stop=stop, num=bins))
            bins.insert(0, 0)
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
        # Save it #
        self.save_plot(fig, axes, **kwargs)
        pyplot.close(fig)
        # For convenience #
        return self

################################################################################
class BinNucleotideDistribution(Graph):
    """Bin total nucleotide size distribution"""
    short_name = 'bins_nucleotide_dist'
    sep        = 'x'
    y_scale    = 'symlog'

    def plot(self, bins=80, **kwargs):
        # Data #
        counts = [sum(map(len, b.contigs)) for b in self.parent.bins]
        # Linear bins in logarithmic space #
        if 'log' in kwargs.get('x_scale', ''):
            start, stop = numpy.log10(1), numpy.log10(max(counts))
            bins = list(numpy.logspace(start=start, stop=stop, num=bins))
            bins.insert(0, 0)
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
    sep        = 'x'
    y_scale    = 'symlog'

    def plot(self, bins=80, **kwargs):
        # Data #
        counts = [b.evaluation.results[''] for b in self.parent.bins]
        # Linear bins in logarithmic space #
        if 'log' in kwargs.get('x_scale', ''):
            start, stop = numpy.log10(1), numpy.log10(max(counts))
            bins = list(numpy.logspace(start=start, stop=stop, num=bins))
            bins.insert(0, 0)
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