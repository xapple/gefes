# Built-in modules #

# Internal modules #
from gefes.graphs import Graph

# Third party modules #
import pandas
from matplotlib import pyplot

# Constants #
__all__ = ['CleanOutcome']

################################################################################
class CleanOutcome(Graph):
    """The results of the cleaning procedure"""
    short_name = 'clean_outcome'

    def plot(self):
        # Data #
        rows = [pool.short_name for pool in reversed(self.parent)]
        columns = ['Kept', 'Discarded']
        data = [(len(p.cleaner), len(p)-len(p.cleaner)) for p in reversed(self.parent.pools)]
        self.frame = pandas.DataFrame(data, index=rows, columns=columns)
        # Plot #
        fig = pyplot.figure()
        axes = self.frame.plot(kind='barh', stacked=True)
        fig = pyplot.gcf()
        # Other #
        axes.set_title('Cleaning result (sickle is used)')
        axes.set_xlabel('Number of paired reads')
        axes.xaxis.grid(True)
        axes.yaxis.grid(False)
        # Save it #
        self.save_plot(fig, axes, sep=('x'), left=0.1, right=0.96)
        self.frame.to_csv(self.csv_path)