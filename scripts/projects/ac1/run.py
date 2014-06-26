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
# Clean the pools #
for p in gefes.projects['acI']: p.clean_reads()
gefes.projects['acI'].graphs[0].plot()

#assemble
gefes.projects['acI'].runner.run_slurm(steps=[{'assemble':{}}], partition='large', machines=64, cores=1024, time='1-00:00:00')

#index
gefes.projects['acI'].assembly.index()
#on a dedicated machine
for p in gefes.projects['acI']: p.runner.run_slurm(steps=[{'map_reads':{}}], time='12:00:00',project='default')

clusty = { 'type' : 'GefesKMeans', kwargs : {'nb' : 6 , 'max_freq' : 0.1 , 'min_length' : 1000 , 'transform' = 'rank'}}
gefes.projects['alinen'].binner.new('test', clusty)
gefes.projects['alinen'].binner['test'].run(GefesKMeans(nb=6,max_freq=0.1,min_length=1000,transform='rank'))
