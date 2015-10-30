#!/usr/bin/env python2
# -*- coding: utf-8 -*-

"""
A script that searches for specific pfams in a specific project.
"""

# Modules #
import os, gefes
from gefes.annotation.hmmer import Hmmer
from seqsearch.databases.pfam import SpecificFamily
from plumbing.tmpstuff import new_temp_path
from plumbing.autopaths import DirectoryPath, AutoPaths
from tqdm import tqdm
from fasta import FASTA

# Constants #
home = os.environ['HOME'] + '/'
base_dir = home + 'test/special_pfams/'

###############################################################################
print "Loading."
proj = gefes.projects['soda_rerun'].load()
samples = proj.samples
for s in samples: s.load()
bins = proj.merged.results.binner.results.bins
faa = FASTA(base_dir + 'all_proteins.faa')

###############################################################################
if not faa.exists:
    print "Regrouping bins."
    temp = FASTA(new_temp_path())
    temp.create()
    for b in tqdm(bins): temp.add(b.faa)
    temp.close()
    temp.move_to(faa)

###############################################################################
class CustomPfamSearch(object):
    """When you are interested in having an HMM 'database' with only
    one specific Pfam in it."""

    all_paths = """
    /model.hmm
    /seq_hits.txt
    """

    def __init__(self, fam_name):
        # Attributes #
        self.fam_name = fam_name
        self.pfam = SpecificFamily(self.fam_name)
        # Directory #
        self.base_dir = DirectoryPath(base_dir + self.fam_name + '/')
        self.p        = AutoPaths(self.base_dir, self.all_paths)

    @property
    def search(self): return Hmmer(faa, self.base_dir, self.pfam.hmm_db)

    @property
    def results(self):
        if not self.search:
            print "Running search."
            self.search.run(cpus=4)
        return self.search.results

    @property
    def fasta(self):
        """The fasta file containing the preficted proteins that recieved
        an annotation."""
        fasta = FASTA(self.p.fasta)
        if not fasta:
            fasta.create()
            for gene in self.filtered_genes: fasta.add_str(str(gene), name=gene.name)
            fasta.close()
        return fasta

    @property
    def alignment(self):
        """The fasta file aligned with muscle and filtered with gblocks"""
        muscle    = AlignedFASTA(self.p.muscle)
        alignment = AlignedFASTA(self.p.aln)
        if not alignment:
            self.fasta.align(muscle)
        return alignment

    @property
    def tree(self):
        """The path to the tree built with raxml"""
        tree = FilePath(self.p.tree_dir + 'RAxML_bestTree.tree')
        if not tree.exists:
            self.alignment.build_tree(new_path    = self.p.tree_dir,
                                      seq_type    = self.analysis.seq_type,
                                      num_threads = self.analysis.num_threads,
                                      free_cores  = 0,
                                      keep_dir    = True)
        return tree

###############################################################################
families = ('PF00151.15', 'PF00150.14', 'PF12876.3', 'PF00128.20')
f = families[0]

searches = [CustomPfamSearch(x) for x in families]
s = searches[0]

###############################################################################
#for s in searches:
#    print "Doing the job for family '%s'" % s.f
#    s.run()
