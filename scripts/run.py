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

###############################################################################
# Just one pool via slurm #
gefes.projects['test']['run000-pool01'].run_slurm()
gefes.projects['humic']['run001-pool01'].run_slurm()

# Just one function for one pool #
gefes.projects['test'][0].runner(steps=[{'clean_reads':{}}], threads=False)

###############################################################################
# Clean the pools #
for p in gefes.projects['test']: p.clean_reads()
for p in gefes.projects['humic']: p.runner.run_slurm(steps=[{'clean_reads':{}}])

# The clean graphs #
for p in gefes.projects['test']: p.make_graphs()
for p in gefes.projects['test']: p.runner.run_slurm(steps=[{'make_plots':{}}])
for p in gefes.projects['humic']: p.runner.run_slurm(steps=[{'make_plots':{}}])

# Assemble #
gefes.projects['test'].assemble()
gefes.projects['test'].runner.run_slurm(steps=[{'assemble':{}}])
gefes.projects['humic'].runner.run_slurm(steps=[{'assemble':{}}], cluster='halvan', cores=64)

# The assembly graphs #
gefes.projects['test'].graphs[0].plot()
gefes.projects['humic'].graphs[0].plot()

# Map #
gefes.projects['test'].assembly.index()
for p in gefes.projects['test']: p.runner.run_slurm(steps=[{'map':{}}])
