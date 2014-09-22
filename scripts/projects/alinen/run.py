#!/usr/bin/env python2

"""
A script to contain the procedure for running the test sample.
"""

# Built-in modules #

# Internal modules #
import gefes

################################ Preprocessing ################################
# Constants #
proj = gefes.projects['alinen'].load()
samples = proj.samples

# Manual #
for s in samples:
    s.load()
    s.pair.fwd.fastqc.run()
    s.pair.rev.fastqc.run()
    s.quality_checker.run()
    s.clean.fwd.fastqc.run()
    s.clean.rev.fastqc.run()
    s.clean.fwd.graphs.length_dist.plot()
    s.clean.rev.graphs.length_dist.plot()
    s.pair.fwd.avg_quality
    s.pair.rev.avg_quality
    s.report.generate()

# By SLURM #
for s in samples: s.load().runner.run_slurm(time='04:00:00')

################################### Assembly ##################################
# On Milou #
proj.runner.run_slurm(steps=['assembly41.run'], time='3-00:00:00', constraint='mem512GB', project="g2014124", job_name="alinen_ray_41")
proj.runner.run_slurm(steps=['assembly51.run'], time='3-00:00:00', constraint='mem512GB', project="g2014124", job_name="alinen_ray_51")
proj.runner.run_slurm(steps=['assembly61.run'], time='3-00:00:00', constraint='mem512GB', project="g2014124", job_name="alinen_ray_61")
proj.runner.run_slurm(steps=['assembly71.run'], time='3-00:00:00', constraint='mem512GB', project="g2014124", job_name="alinen_ray_71")

# On Sisu #
proj.runner.run_slurm(steps=['assembly71.run'], machines=64, cores=64*24, time='12:00:00', partition='large', job_name="alinen_ray_71", email=False)

# On Halvan #
proj.runner.run_slurm(steps=['assembly41.run'], time='10-00:00:00', project="b2011035", job_name="alinen_proj_41", cluster='halvan', partition='halvan', cores=64)

################################### Aggregate ##################################
hypo = gefes.groups.favorites.alinen_hypo.load()
meta = gefes.groups.favorites.alinen_meta.load()
epi = gefes.groups.favorites.alinen_epi.load()

meta.runner.run_slurm(steps=[{'assembly41.run':{'threads':False}}], time='7-00:00:00', constraint='mem512GB', project="b2011105", job_name="alinen_epi_41")
epi.runner.run_slurm(steps=[{'assembly41.run':{'threads':False}}], time='7-00:00:00', constraint='mem512GB', project="b2011105", job_name="alinen_hypo_41")
hypo.runner.run_slurm(steps=[{'assembly41.run':{'threads':False}}], time='7-00:00:00', constraint='mem512GB', project="b2011105", job_name="alinen_hypo_41", machines=4, cores=4*16)