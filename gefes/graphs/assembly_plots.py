# Futures #
from __future__ import division

# Built-in modules #
# Internal modules #
from gefes.graphs import Graph
from gefes.groups.pools import Pool
# Third party modules #
import pandas, matplotlib
from matplotlib import pyplot
import scipy.stats
import math

# Constants #
__all__ = ['ContigDist', 'ContigCumsum','ContigDistribution']

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


################################################################################
class ContigDistribution(Graph):
    """Cumulative sum of nuclotides by contig lengths (c.f. brouillon)"""
    short_name = 'contig_distrib'

    def plot(self, pool_list = None,max_freq=None,min_len=None, col_mode="gc_content"):
        # Data #
        
        if pool_list == None : pools = [p.id_name for p in self.parent.aggregate]
        else :
            pools = [p.id_name for p in self.parent.aggregate][pool_list]
        print(pools)
            
        gframe=self.parent.filtered_frame(max_freq=max_freq,min_len=min_len)
        if isinstance(pools,list):
                yvalues = gframe[pools].mean(1)
        else:
                yvalues = gframe[pools]
        title = 'Distribution of Contigs for ' + "/".join(pools)

        xvalues = gframe["length"]
        yvalues = [y+1 for y in yvalues]

        if col_mode is "gc_content":
            zvalues = scipy.stats.rankdata(gframe["gc_content"])
            quantiles = [int(math.floor(i*(len(zvalues)-1.0)/5.0))for i in range(0,6)]
            breaks= [(j,i-1) for i,j in zip(zvalues,gframe["gc_content"]) if i-1 in quantiles]
            if len(breaks) != len(quantiles):
                    quantiles = list(set([b[1] for b in breaks if int(b[1]) in quantiles] ))
                    quantiles.sort()
            breaks=list(set([b[0] for b in breaks]))
            breaks.sort()
            breaks=["%.2f" % b for b in breaks]

        else:
#            if not self.parent.aggregate.binner.clusterer.kmeans.coverage_clusters: self.parent.aggregate.binner.clusterer.kmeans.run()
            zvalues = self.parent.aggregate.binner.clusterer.kmeans.tetras_clusters

        

        
        # Plot #
        cm = pyplot.cm.get_cmap('gist_rainbow')
        fig = pyplot.figure()

        pyplot.scatter(xvalues, yvalues, c = zvalues, cmap = cm)
        
        if col_mode is "gc_content":
            cax = pyplot.cm.ScalarMappable(cmap=cm)
            cax.set_array(zvalues)
            cbar = fig.colorbar(cax,ticks=quantiles)
            cbar.ax.set_yticklabels(breaks)


        axes = pyplot.gca()
        axes.set_title(title)
        axes.set_xlabel('Length of contigs in nucleotides')
        axes.set_ylabel('Coverage')
        axes.set_xscale('log')
        axes.set_yscale('log')
        axes.yaxis.grid(True)
        axes.xaxis.grid(True)
        
        # Save it #
        self.save_plot(fig, axes, left=0.09, height=12.0)
        pyplot.close(fig)





