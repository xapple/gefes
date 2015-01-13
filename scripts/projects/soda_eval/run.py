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

# Three which fail with memory errors #
fail = [proj['ss10'], proj['ss12'], proj['zl01'], proj['zl12']]
for s in fail: s.runner.run_slurm(steps=['mono_mapper.run'], machines=1, cores=32, time='1-00:00:00', partition='hugemem')

################################ Solo-Annotation ###############################
from tqdm import tqdm
all_contigs = [c for s in proj.samples for c in s.contigs]
for c in tqdm(all_contigs): c.annotation.run()
