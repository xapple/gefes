#!/usr/bin/env python

"""
A script to contain examples commands
for running the pipeline.
"""

# Don't run it #
import sys
sys.exit("Copy paste the commands you want in ipython, don't run this script.")

# Modules #
import gefes

# Just one pool via slurm #
gefes.projects['test']['run000-pool01'].run_slurm()
gefes.projects['humic']['run001-pool01'].run_slurm()

# Just one function for one pool #
pj = gefes.projects['test']; p = pj[0]; p(steps=[{'clean_reads':{}}], threads=False)

# Clean the pools #
for p in gefes.projects['test']: p.cleaner.clean()
[p.runner.run_slurm(steps=[{'clean_reads':{}}]) for p in gefes.projects['test']]
# The clean graphs #
gefes.projects['test'].graphs[0].plot()
for p in gefes.projects['test']: p.graphs[0].plot()

# Assemble #
gefes.projects['test'].assemble()
gefes.projects['test'].runner.run_slurm(steps=[{'assemble':{}}])
# The assembly graphs #

# Map #
for p in gefes.projects['test']: p.mapper.map()

