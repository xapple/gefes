#!/usr/bin/env python2

"""
A script to contain the procedure for running the alinen project.
"""

# Built-in modules #

# Internal modules #
import gefes

# Third party modules #
from tqdm import tqdm

#################################### Load #####################################
# One project #
proj = gefes.projects['alinen'].load()
samples = proj.samples
for s in samples: s.load()

# Don't run it #
import sys
sys.exit("Copy paste the commands you want in ipython, don't run this script.")

################################ Status report ################################
# How far did we run things #
for s in samples: print "Raw:",                s, bool(s.pair)
for s in samples: print "First QC:",           s, bool(s.pair.fwd.fastqc)
for s in samples: print "Cleaned:",            s, bool(s.quality_checker.results)
for s in samples: print "Second QC:",          s, bool(s.clean.fwd.fastqc)
for s in samples: print "Initial taxa:",       s, bool(s.kraken)
for s in samples: print "Mono-assembly:",      s, bool(s.assembly)
for k,v in proj.assemblies.items(): print "Co-assembly %i:"%k, proj, bool(v)
for s in samples: print "Mono-mapping:",       s, bool(s.mono_mapper.p.coverage)
print                   "Merged assembly:",  proj, bool(proj.merged.results)
for s,a,m in ((s,a,m) for a,m in s.mappers.items() for s in samples): print "Map %s to %s:"%(s,a), bool(m.p.coverage)
for k,v in proj.assemblies.items(): print "Binning %i:"%k, proj, bool(v.results.binner.results)
print                   "Merged binning:", bool(proj.merged.results.binner.results)

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
for s in proj.samples:
    s.runner.run_slurm(steps=['assembly.run'], machines=3, cores=3*24, time='12:00:00', partition='small')

################################ Solo-Mapping #################################
for s in samples: print s, s.mono_mapper.run()

############################### Solo-Prokka #################################
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

################################# Co-Mapping ##################################
for s in samples: s.mapper.run()

################################# Concoct ##################################
proj.assembly.results.binner.run()

################################# Phylosift ##################################
for c in proj.assembly.results.contigs: c.taxonomy.run()

############################## Different kmers ##################################
# Run assembly #
proj.runner.run_slurm(steps=['assembly_51.run'], machines=42, cores=42*24,
                      time='36:00:00', partition='large', job_name="alinen_ray_51")
proj.runner.run_slurm(steps=['assembly_61.run'], machines=42, cores=42*24,
                      time='36:00:00', partition='large', job_name="alinen_ray_61")
proj.runner.run_slurm(steps=['assembly_81.run'], machines=42, cores=42*24,
                      time='36:00:00', partition='large', job_name="alinen_ray_81")

# Merge with Newbler #
params = dict(machines=1, cores=24, time='14-00:00:00', partition='longrun', constraint='hsw', memory=120000)
proj.runner.run_slurm(steps=['merged.run'], job_name="gefes_newbler", **params)

# Run mapping #
params = dict(machines=1, cores=1, threads=24, time='14-00:00:00', partition='longrun', constraint='hsw', memory=120000)
for s in samples: s.runner.run_slurm(steps=['mapper_51.run'], job_name=s.name + "_co_51_map", **params)
for s in samples: s.runner.run_slurm(steps=['mapper_61.run'], job_name=s.name + "_co_61_map", **params)
for s in samples: s.runner.run_slurm(steps=['mapper_71.run'], job_name=s.name + "_co_71_map", **params)
for s in samples: s.runner.run_slurm(steps=['mapper_81.run'], job_name=s.name + "_co_81_map", **params)
for s in samples: s.runner.run_slurm(steps=['mapper_merged.run'], job_name=s.name + "_merge_map", **params)

# Run Concoct #
params = dict(machines=1, cores=1, time='14-00:00:00', partition='longrun', constraint='hsw', memory=120000)
proj.runner.run_slurm(steps=['assembly_51.results.binner.run'], job_name="concot_51", **params)
proj.runner.run_slurm(steps=['assembly_61.results.binner.run'], job_name="concot_61", **params)
proj.runner.run_slurm(steps=['assembly_71.results.binner.run'], job_name="concot_71", **params)
proj.runner.run_slurm(steps=['assembly_81.results.binner.run'], job_name="concot_81", **params)
proj.runner.run_slurm(steps=['merged.results.binner.run'], job_name="concot_merged", **params)

# Run CheckM #
params = dict(machines=1, cores=1, time='3-00:00:00', partition='serial', constraint='hsw', memory=124000)
proj.runner.run_slurm(steps=['assembly_51.results.binner.results.run_all_bin_eval'], job_name="checkm_51", **params)
proj.runner.run_slurm(steps=['assembly_61.results.binner.results.run_all_bin_eval'], job_name="checkm_61", **params)
params = dict(cores=10, time='24:00:00', partition='hugemem', memory=400000)
proj.runner.run_slurm(steps=['assembly_71.results.binner.results.run_all_bin_eval'], job_name="checkm_71", **params)
proj.runner.run_slurm(steps=['assembly_81.results.binner.results.run_all_bin_eval'], job_name="checkm_81", **params)
params = dict(machines=1, cores=1, memory=124000, time='1-00:00:00', partition='serial', constraint='hsw')
proj.runner.run_slurm(steps=['merged.results.binner.results.run_all_bin_eval'], job_name="checkm_merged", **params)
