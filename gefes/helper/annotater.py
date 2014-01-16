# Futures #
from __future__ import division

# Built-in modules #

# Internal modules #
from gefes.common import flatten
from gefes.common.autopaths import AutoPaths
from gefes.common.cache import property_cached

from gefes.helper.contig import Contig
from gefes.fasta.single import FASTA
from gefes.common.slurm import nr_threads

import sh

class Annotation(object):
    """The genomic annotation from one annotation tool"""
    
    all_paths = """
/annotation.gff
/short_form.txt
/Phyl_AMPHORA/
/prodigal_blastn/
"""
    tools=['prodigal:blastn','Phyl_AMPHORA']

    def __init__(self,contig,tool,parameters):
        # Save parent #
        self.parent, self.contig = contig, contig
        self.base_dir = self.parent.base_dir
        self.p = AutoPaths(self.base_dir, self.all_paths)
        if tool in tools: self.tool=tool
        else: self.tool=tools[0]
        sel.parameters=parameters


    def annotate_prodigal_blastn(self):
        # folder where to output intermediate files #
        out_dir = self.p.output_dir
        command = sh.Command('prodigal')
        stats = command('-k', 81, '-o', out_dir, *pairs)
        # Print the report #
        with open(self.p.report, 'w') as handle: handle.write(str(stats))

        
