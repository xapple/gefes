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
p = gefes.projects['strain_mock'].load()
print p.assembly.results.contigs_fasta.count