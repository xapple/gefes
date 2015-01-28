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
# How far did we run things #
for s in samples: print "Raw:",                s, bool(s.pair)
for s in samples: print "First QC:",           s, bool(s.pair.fwd.fastqc.results)
for s in samples: print "Cleaned:",            s, bool(s.quality_checker.results)
for s in samples: print "Second QC:",          s, bool(s.clean.fwd.fastqc.results)
for s in samples: print "Initial taxa:",       s, bool(s.kraken.results)
for s in samples: print "Solo-assembly:",      s, bool(s.assembly.results)
for s in samples: print "Mono-mapping:",       s, bool(s.mono_mapper.results)
for s in samples: print "Map to co-assebmly:", s, bool(s.mapper.results)

# Project #

# Logs #
for s in samples: print "Logs:",               s, list(s.p.logs_dir.contents)
print "Logs:", proj, list(proj.p.logs_dir.contents)

# Report #
for s in samples: s.report.generate()

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

################################### Kraken ####################################
for s in samples: print s, s.kraken.run()

############################### Solo-Assembly #################################
for s in proj.samples: s.runner.run_slurm(steps=['assembly.run'], machines=3, cores=3*24, time='12:00:00', partition='small')

################################ Solo-Mapping #################################
for s in samples: print s, s.mono_mapper.run()

############################### Solo-Prokka #################################
from tqdm import tqdm
all_contigs = [c for s in proj.samples for c in s.contigs]
for c in tqdm(all_contigs): c.annotation.run()

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