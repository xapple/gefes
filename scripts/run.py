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

# Just one function for one pool #
pj = gefes.projects['test']; p = pj[0]; p(steps=[{'assemble':{}}], threads=False)

# Just one pool via slurm #
gefes.projects['test']['run000-pool01'].run_slurm()
gefes.projects['humic']['run001-pool01'].run_slurm()