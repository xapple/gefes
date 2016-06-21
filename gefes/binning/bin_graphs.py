# Built-in modules #

# Internal modules #
from plumbing.graphs import Graph

# Third party modules #
import numpy, scipy
from matplotlib import pyplot
from Bio.SeqUtils import GC

# Constants #
__all__ = ['GcFractionOnTotalCoverage']

################################################################################
class GcFractionOnTotalCoverage(Graph):
    """Bin number of contigs distribution"""
    short_name = 'gc_x_totcov'
    y_grid     = True
    title      = 'Contig GC content against total coverage'
    x_label    = 'GC content'
    y_label    = 'Total coverage'
    y_scale    = 'symlog'

    def plot(self, **kwargs):
        # Data #
        x = [GC(c.record.seq) for c in self.parent.contigs]
        y = [self.parent.binner.coverage_matrix.loc[c.name].sum() for c in self.parent.contigs]
        # Plot #
        fig = pyplot.figure()
        pyplot.plot(x, y, 'go')
        axes = pyplot.gca()
        # Save it #
        self.save_plot(fig, axes, **kwargs)
        pyplot.close(fig)
        # For convenience #
        return self
