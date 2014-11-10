#!/usr/bin/env python2

"""
A script to contain examples commands for running the pipeline.
"""

# Don't run it #
import sys
sys.exit("Copy paste the commands you want in ipython, don't run this script.")

# Modules #
import gefes

###############################################################################
from plumbing import slurm
slurm.nr_threads = 32

p = gefes.projects['strain_mock'].load()
s = p.first.load()

print p.assembly.results.contigs_fasta.count
s.mapper.run()
print s.mapper.results.coverage