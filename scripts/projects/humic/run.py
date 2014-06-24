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
# Just one pool via slurm #
gefes.projects['humic']['run001-pool01'].run_slurm()

# Just one function for one pool #
gefes.projects['humic'][0].runner(steps=[{'clean_reads':{}}], threads=False)

###############################################################################
# Fastqc on the pools #
for p in gefes.projects['humic']: p.pair.fastqc(directory=p.p.fastqc_dir)
# Clean the pools #
for p in gefes.projects['humic']: p.runner.run_slurm(steps=[{'clean_reads':{}}])
# The clean graphs #
for p in gefes.projects['humic']: p.runner.run_slurm(steps=[{'make_plots':{}}])
gefes.projects['humic'].graphs[0].plot()
# Fastqc on the result #
for p in gefes.projects['humic']: p.cleaner.pair.fastqc(directory=p.cleaner.p.fastqc_dir)


# Assemble on kalykl #
gefes.projects['test'].runner.run_slurm(steps=[{'assemble':{}}])
# Assemble on halvan #
gefes.projects['humic'].runner.run_slurm(steps=[{'assemble':{}}], cluster='halvan', cores=64, time='6-12:00:00')
# Assemble on sisu #
gefes.projects['humic'].runner.run_slurm(steps=[{'assemble':{}}], partition='large', machines=64, cores=1024, time='1-00:00:00')
# The assembly graphs #
gefes.projects['humic'].assembly.graphs[0].plot()

# Index the result #
gefes.projects['humic'].assembly.index()
# Map the reads #
for p in gefes.projects['humic']: p.runner.run_slurm(steps=[{'map_reads':{}}], time='12:00:00')

# Binning frame #
gefes.projects['humic'].binner.export_frame()
# Clustering #
gefes.projects['humic'].binner.clusterer.run()