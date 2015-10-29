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
from plumbing.autopaths import DirectoryPath, FilePath

# Constants #
home = os.environ['HOME'] + '/'

###############################################################################
print "Loading"
proj = gefes.projects['soda_rerun'].load()
samples = proj.samples
for s in samples: s.load()
good_contigs = proj.merged.results.binner.results.good_contigs

print "Extracting families"
small_db = FilePath(home + "databases/pfam/temp_small.hmm")
families = ('PF00151.15', 'PF00150.14', 'PF12876.3', 'PF00128.20')
fam_file = TmpFile.from_string('\n'.join(families))
sh.hmmfetch('-o', small_db, '-f', pfam.hmm_db, fam_file)
assert small_db

#print "Directories"
#dirs     = [DirectoryPath(c.base_dir + 'special_pfams/')                               for c in good_contigs]
#for d in dirs: d.create(safe=True)
#
#print "Searching"
#searches = [Hmmer(c.proteins.results.faa, c.base_dir + 'special_pfams/', small_db) for c in good_contigs]
#for s in searches: pass#