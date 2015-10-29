#!/usr/bin/env python2
# -*- coding: utf-8 -*-

"""
A script that searches for specific pfams in a specific project.
"""

# Modules #
import os, sh, gefes
from gefes.annotation.hmmer import Hmmer
from seqsearch.databases.pfam import pfam
from plumbing.tmpstuff import new_temp_path
from plumbing.autopaths import DirectoryPath, FilePath
from tqdm import tqdm
from fasta import FASTA

# Constants #
home = os.environ['HOME'] + '/'

###############################################################################
print "Loading."
proj = gefes.projects['soda_rerun'].load()
samples = proj.samples
for s in samples: s.load()

###############################################################################
print "Regrouping bins."
base_dir    = home + 'test/special_pfams/'
bins        = proj.merged.results.binner.results.bins
faa         = FASTA(base_dir + 'all_proteins.faa')
if not faa.exists:
    temp = FASTA(new_temp_path())
    temp.create()
    for b in tqdm(bins): temp.add(b.faa)
    temp.close()
    temp.move_to(faa)

###############################################################################
class CustomPfamSearch(object):
    def __init__(self, f):
        self.f = f
        self.dir = DirectoryPath(base_dir + f + '/')
        self.dir.create(safe=True)

    @property
    def hmm_db(self):
        hmm_db = FilePath(self.dir + 'db.hmm')
        if not hmm_db.exists:
            print sh.hmmfetch('-o', hmm_db, pfam.hmm_db, f)
            assert hmm_db
        return hmm_db

    @property
    def search(self): return Hmmer(faa, self.dir, self.hmm_db)

    @property
    def search_results(self):
        if not self.search:
            print "Running search."
            self.search.run(cpus=4)
        return self.search.results

    @property
    def tree(self): return Hmmer(faa, self.dir, self.hmm_db)

###############################################################################
families = ('PF00151.15', 'PF00150.14', 'PF12876.3', 'PF00128.20')
f = families[0]

searches = [CustomPfamSearch(x) for x in families]
s = searches[0]

###############################################################################
#for s in searches:
#    print "Doing the job for family '%s'" % s.f
#    s.run()
