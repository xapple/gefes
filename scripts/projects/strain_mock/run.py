#!/usr/bin/env python2

"""
A script to contain examples commands for running the pipeline.
"""

# Don't run it #
import sys
sys.exit("Copy paste the commands you want in ipython, don't run this script.")

# Modules #
import gefes
from fasta import FASTA

###############################################################################
# Custom threads #
from plumbing import slurm
slurm.nr_threads = 32

# Project and sample objects #
p = gefes.projects['strain_mock'].load()
for s in p: s.load()

# Testing the pipeline #
print p.assembly.results.contigs_fasta.count
s.mapper.run()
print s.mapper.results.coverage

# Filter contigs #
FASTA(p.assembly.p.Contigs).extract_length(new_path=p.assembly.p.filtered, lower_bound=1000)

# Compute all coverages #
for sample in p: print """python -c "import gefes; p = gefes.projects['strain_mock'].load(); i=%s; s = p[i].load(); s.mapper.results.statistics" &""" % sample.num

# Run concoct #
p.binner.run()
b = p.binner.results.bins[0]

# Annotate #
b.annotation.run()