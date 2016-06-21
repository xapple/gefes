# Built-in modules #

# Internal modules #
from plumbing.graphs import Graph

# Third party modules #
import numpy, scipy
from matplotlib import pyplot

# Constants #
__all__ = ['BinContigDistribution', 'BinNucleotideDistribution', 'BinGenesOnPfams',
           'BinBasesOnGenes']

################################################################################
class BinContigDistribution(Graph):
    """Bin number of contigs distribution"""
    short_name = 'bins_contig_dist'
    sep        = 'y'
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
    sep        = 'y'
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
class BinGenesOnPfams(Graph):
    """bin_genes_x_pfams"""

    short_name = 'bin_genes_x_pfams'
    title      = 'Bins: predicted number of genes VS the number of predicted PFAMs'
    x_label    = 'Number of predicted genes in bin'
    y_label    = 'Number of distinct predicted PFAMs in bin'
    sep        = 'x'
    y_grid     = True
    left       = 0.1

    def plot(self, **kwargs):
        # Data #
        x      = [b.faa.count                 for b in self.parent.bad_bins]
        y      = [len(b.pfams.distinct_pfams) for b in self.parent.bad_bins]
        x_good = [b.faa.count                 for b in self.parent.good_bins]
        y_good = [len(b.pfams.distinct_pfams) for b in self.parent.good_bins]
        # Plot #
        fig = pyplot.figure()
        pyplot.plot(x,      y,      'bo', label='Other bins')
        pyplot.plot(x_good, y_good, 'ro', label='Bins marked good')
        # Regression #
        slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(x+x_good,y+y_good)
        regress_y = intercept + slope * numpy.array(x)
        pyplot.plot(x, regress_y, 'k-', label='Linear regression R=%.2g' % r_value)
        pyplot.legend()
        axes = pyplot.gca()
        # Save it #
        self.save_plot(fig, axes, **kwargs)
        pyplot.close(fig)
        # For convenience #
        return self

################################################################################
class BinBasesOnGenes(Graph):
    """bin_bps_x_genes"""

    short_name = 'bin_bps_x_genes'
    title      = 'Bins: length in base pairs VS the number of predicted genes'
    x_label    = 'Cumulative length in base pairs of all contigs in bin'
    y_label    = 'Number of predicted genes in bin'
    sep        = 'x'
    y_grid     = True
    left       = 0.1

    def plot(self, **kwargs):
        # Data #
        x      = [sum(map(len, b.contigs)) for b in self.parent.bad_bins]
        y      = [b.faa.count              for b in self.parent.bad_bins]
        x_good = [sum(map(len, b.contigs)) for b in self.parent.good_bins]
        y_good = [b.faa.count              for b in self.parent.good_bins]
        # Plot #
        fig = pyplot.figure()
        pyplot.plot(x,      y,      'bo', label='Other bins')
        pyplot.plot(x_good, y_good, 'ro',  label='Bins marked good')
        # Regression #
        slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(x+x_good,y+y_good)
        regress_y = intercept + slope * numpy.array(x)
        pyplot.plot(x, regress_y, 'k-', label='Linear regression R=%.2g' % r_value)
        pyplot.legend()
        axes = pyplot.gca()
        # Save it #
        self.save_plot(fig, axes, **kwargs)
        pyplot.close(fig)
        # For convenience #
        return self
