#!/usr/bin/env python2

"""
A script to contain the procedure for running the test sample.
"""

# Built-in modules #

# Internal modules #
import gefes

# Constants #
proj = gefes.projects['alinen'].load()
samples = proj.samples
for s in samples: s.load()

################################ Status report ################################
for s in samples: print "Raw:",              s, bool(s.pair)
for s in samples: print "First QC:",         s, bool(s.pair.fwd.fastqc.results)
for s in samples: print "Cleaned:",          s, bool(s.clean)
for s in samples: print "Quality checked:",  s, bool(s.quality_checker.results)
for s in samples: print "Initial taxa fwd:", s, bool(s.clean.fwd.kraken.results)
for s in samples: print "Initial taxa rev:", s, bool(s.clean.rev.kraken.results)
for s in samples: print "Solo-assebmly:",    s, bool(s.assembly.results)
for s in samples: print "Mono-mapping:",     s, bool(s.mono_mapper.results)

for s in samples: print "Map to co-assebmly:", s, bool(s.clean.rev.kraken.results)

################################ Preprocessing ################################
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
for s in samples: s.runner.run_slurm(time='04:00:00')

################################# Co-Assembly #################################
# On Milou #
proj.runner.run_slurm(steps=['assembly.run'], time='3-00:00:00', constraint='mem512GB', project="g2014124",
                      job_name="alinen_ray")

# On Sisu #
proj.runner.run_slurm(steps=['assembly.run'], machines=42, cores=42*24, time='36:00:00', partition='large',
                      job_name="alinen_ray", email=False)

# On Halvan #
proj.runner.run_slurm(steps=['assembly.run'], time='10-00:00:00', project="b2011035",
                      job_name="alinen_ray", cluster='halvan', partition='halvan', cores=64)

################################# Aggregates ##################################
hypo = gefes.groups.favorites.alinen_hypo.load()
meta = gefes.groups.favorites.alinen_meta.load()
epi  = gefes.groups.favorites.alinen_epi.load()

meta.runner.run_slurm(steps=[{'assembly41.run':{'threads':False}}], time='7-00:00:00', constraint='mem512GB', project="b2011105", job_name="alinen_epi_41")
epi.runner.run_slurm(steps=[{'assembly41.run':{'threads':False}}], time='7-00:00:00', constraint='mem512GB', project="b2011105", job_name="alinen_hypo_41")
hypo.runner.run_slurm(steps=[{'assembly41.run':{'threads':False}}], time='7-00:00:00', constraint='mem512GB', project="b2011105", job_name="alinen_hypo_41", machines=4, cores=4*16)