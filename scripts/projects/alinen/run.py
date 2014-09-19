#!/usr/bin/env python2

"""
A script to contain the procedure for running the test sample.
"""

# Built-in modules #

# Internal modules #
import gefes

################################ Preprocessing ################################
proj = gefes.projects['alinen']
samples = proj.samples
for s in samples:
    s.load()
    s.pair.fwd.fastqc.run()
    s.pair.rev.fastqc.run()
    s.quality_checker.run()
    s.clean.fwd.fastqc.run()
    s.clean.rev.fastqc.run()
    s.clean.fwd.graphs['LengthDist'].plot()
    s.clean.rev.graphs['LengthDist'].plot()
    s.pair.fwd.avg_quality
    s.pair.rev.avg_quality
    s.report.generate()

for s in samples: s.runner.run_slurm(partition='test_large', time='04:00:00')

################################### Assembly ##################################
proj = gefes.projects['alinen'].load()

proj.runner.run_slurm(steps=[{'assembly41.run':{'threads':False}}], time='3-00:00:00', constraint='mem512GB', project="g2014124", job_name="alinen_ray_41")
proj.runner.run_slurm(steps=[{'assembly51.run':{'threads':False}}], time='3-00:00:00', constraint='mem512GB', project="g2014124", job_name="alinen_ray_51")
proj.runner.run_slurm(steps=[{'assembly61.run':{'threads':False}}], time='3-00:00:00', constraint='mem512GB', project="g2014124", job_name="alinen_ray_61")
proj.runner.run_slurm(steps=[{'assembly71.run':{'threads':False}}], time='3-00:00:00', constraint='mem512GB', project="g2014124", job_name="alinen_ray_71")

proj.runner.run_slurm(steps=[{'assembly41.run':{'threads':False}}], time='10-00:00:00', project="b2011035", job_name="alinen_proj_41", cluster='halvan', partition='halvan', cores=64)

################################### Aggregate ##################################
hypo = gefes.groups.favorites.alinen_hypo.load()
meta = gefes.groups.favorites.alinen_meta.load()
epi = gefes.groups.favorites.alinen_epi.load()

meta.runner.run_slurm(steps=[{'assembly41.run':{'threads':False}}], time='7-00:00:00', constraint='mem512GB', project="b2011105", job_name="alinen_epi_41")
epi.runner.run_slurm(steps=[{'assembly41.run':{'threads':False}}], time='7-00:00:00', constraint='mem512GB', project="b2011105", job_name="alinen_hypo_41")
hypo.runner.run_slurm(steps=[{'assembly41.run':{'threads':False}}], time='7-00:00:00', constraint='mem512GB', project="b2011105", job_name="alinen_hypo_41", machines=4, cores=4*16)