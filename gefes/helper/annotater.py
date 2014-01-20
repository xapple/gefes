# Futures #
from __future__ import division

# Built-in modules #

# Internal modules #
from gefes.common import flatten
from gefes.common.autopaths import AutoPaths
from gefes.common.cache import property_cached

from gefes.fasta.single import FASTA
from gefes.common.slurm import nr_threads

import sh

class Annotation(object):
    """The genomic annotation from one annotation tool"""
    
    all_paths = """
/annotation.gff
/short_form.txt
/Phyl_AMPHORA/
"""
    tools=['Phyl_AMPHORA']

    def __init__(self,contig,tool,parameters):
        # Save parent #
        self.parent, self.contig = contig, contig
        self.base_dir = self.parent.base_dir
        self.p = AutoPaths(self.base_dir, self.all_paths)
#        if tool in tools: self.tool=tool
#        else: self.tool=tools[0]
        self.parameters=parameters

###############################################################################
