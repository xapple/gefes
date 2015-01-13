#!/usr/bin/env python2

"""
A script to contain the procedure for running the soda evaluation project.
"""

# Built-in modules #

# Internal modules #
import gefes

################################ Preprocessing ################################
# Constants #
proj = gefes.projects['soda_eval'].load()
samples = proj.samples

# Clean #
for s in proj.samples: s.load()
for s in proj.samples: s.runner.run_slurm()

################################### Solo-Assembly ##################################
for s in proj.samples: s.runner.run_slurm(steps=['assembly.run'], machines=3, cores=3*24, time='12:00:00', partition='small')

################################### Solo-Mapping ##################################
for s in proj.samples: s.runner.run_slurm(steps=['mono_mapper.run'], machines=1, cores=16, time='3-00:00:00', partition='serial')
