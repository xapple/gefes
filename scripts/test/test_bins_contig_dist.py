#!/usr/bin/env python2

"""
A script to test the bins_contig_dist graph.
"""

# Copy #
import sh
source = "warwick-node:/home/lucas/GEFES/reports/"
dest   = "~/Desktop/contigs.fasta"
sh.rsync('-avz', source, dest)

# Do it #
from gefes.binning.graphs import BinContigDistribution

class Bin(object): pass
bin = Bin()
for graph in graphs.__all__:
    cls = getattr(graphs, graph)
    setattr(result, cls.short_name, cls(self))


graph = BinContigDistribution()
