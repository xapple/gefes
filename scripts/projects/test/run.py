#!/usr/bin/env python2

"""
A script to contain the procedure for running the test sample.
"""

# Built-in modules #

# Internal modules #
import gefes

################################ Preprocessing ################################
# Create #
execfile("~/repos/gefes/scripts/projects/alinen/generate_test.py")

# Constants #
proj = gefes.projects['test'].load()
samples = proj.samples

# Clean #
proj.run_samples()
for s in samples: s.runner.run_slurm(partition='test', time='00:15:00')

################################### Solo-Assembly ##################################
sample = proj[0].load()
sample.runner.run_slurm(steps=['assembly.run'], machines=3, cores=3*24, time='12:00:00', partition='small', job_name="test1_ray_71")


################################### Co-Assembly ##################################
# On Milou #
proj.runner.run_slurm(steps=['assembly.run'], time='00:15:00', qos='short', partition='devel', job_name="test_ray_41")

# On Sisu #
proj.runner.run_slurm(steps=['assembly.run'], machines=3, cores=3*24, time='00:30:00', partition='test', job_name="test_ray_41")
proj.runner.run_slurm(steps=['assembly.run'], machines=3, cores=3*24, time='04:00:00', partition='test-large', job_name="test_ray_41")

# On Halvan #
proj.runner.run_slurm(steps=['assembly41.run'], cores=16, time='00:15:00', project="b2011035", job_name="test_41", cluster='halvan', partition='halvan')

################################### Aggregate ##################################
a = gefes.groups.favorites.test_agg
a.load()
a.runner.run_slurm(steps=[{'assembly.run':{}}], time='00:15:00', partition='devel', threads=False)

################################### Mapping ##################################
sample = proj[0].load()

# On taito #
sample.runner.run_slurm(steps=['mono_mapper.run'], machines=1, cores=16, time='00:30:00', partition='test', job_name="test1_bowtie")
