#!/usr/bin/env python2
# -*- coding: utf-8 -*-

"""
A script that searches for specific pfams in a specific project.
"""

# Modules #
import os, sh, gefes
from gefes.annotation.hmmer import Hmmer
from seqsearch.databases.pfam import pfam
from plumbing.tmpstuff import TmpFile

# Constants #
home = os.environ['HOME'] + '/'

###############################################################################
proj = gefes.projects['soda_rerun'].load()
samples = proj.samples
for s in samples: s.load()
good_contigs = proj.merged.results.binner.results.good_contigs

small_db = home + "databases/pfam/temp_small.hmm"
families = ('PF00151', 'PF00150', 'PF12876', 'PF00128')
fam_file = TmpFile.from_string('\n'.join(families))
sh.hmmfetch('-o', small_db, '-f', pfam.hmm_db, fam_file)

#searches = [Hmmer(c.proteins.results.faa, small_db, c.base_dir + 'special_pfams') for c in good_contigs]
