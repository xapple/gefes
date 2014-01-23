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
for p in gefes.projects['test']: p.make_plots()
for p in gefes.projects['humic']: p.runner.run_slurm(steps=[{'make_plots':{}}])
gefes.projects['test'].graphs[0].plot()
gefes.projects['humic'].graphs[0].plot()

# Assemble locally #
gefes.projects['test'].assemble()
# Assemble on kalykl #
gefes.projects['test'].runner.run_slurm(steps=[{'assemble':{}}])
# Assemble on halvan #
gefes.projects['test'].runner.run_slurm(steps=[{'assemble':{}}], cluster='halvan', cores=16, time='15:00')
gefes.projects['humic'].runner.run_slurm(steps=[{'assemble':{}}], cluster='halvan', cores=64, time='6-12:00:00')
# Assemble on sisu #
gefes.projects['test'].runner.run_slurm(steps=[{'assemble':{}}], partition='test', machines=16, cores=256)
gefes.projects['humic'].runner.run_slurm(steps=[{'assemble':{}}], partition='large', machines=64, cores=1024, time='1-00:00:00')
# The assembly graphs #
gefes.projects['test'].assembly.graphs[0].plot()

# Index the result #
gefes.projects['test'].assembly.index()
gefes.projects['humic'].assembly.index()
# Map the reads #
for p in gefes.projects['test']: p.mapper.map()
for p in gefes.projects['test']: p.runner.run_slurm(steps=[{'map_reads':{}}])
for p in gefes.projects['humic']: p.map()
for p in gefes.projects['humic']: p.runner.run_slurm(steps=[{'map_reads':{}}], time='12:00:00')

# Binning frame #
gefes.projects['test'].binner.export_frame()
gefes.projects['humic'].binner.export_frame()
# Clustering #
gefes.projects['test'].binner.clusterer.run()


###############################################################################

# Clean the pools #
for p in gefes.projects['acI']: p.clean_reads()
gefes.projects['acI'].graphs[0].plot()

#assemble
gefes.projects['acI'].runner.run_slurm(steps=[{'assemble':{}}], partition='large', machines=64, cores=1024, time='1-00:00:00')

#index
gefes.projects['acI'].assembly.index()
#on a dedicated machine
for p in gefes.projects['acI']: p.runner.run_slurm(steps=[{'map_reads':{}}], time='12:00:00',project='default')
single_acIs=[gefes.Project('acI'+p.id_name,p,'/homeappl/home/buck/GEFES/views/projects/') for p in gefes.projects['acI']]
for p in gefes.projects['acI']: p.runner.run_slurm(steps=[{'assemble':{}}], partition='large', machines=64, cores=1024, time='1-00:00:00')
