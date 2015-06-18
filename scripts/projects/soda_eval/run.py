#!/usr/bin/env python2

"""
A script to contain the procedure for running the soda evaluation project.
"""

# Built-in modules #

# Internal modules #
import gefes

#################################### Load #####################################
# One project #
proj = gefes.projects['soda_eval'].load()
samples = proj.samples
for s in samples: s.load()

# Don't run it #
import sys
sys.exit("Copy paste the commands you want in ipython, don't run this script.")

################################## Meta-data ##################################
# Print number of sequences #
for s in samples: print s.pair.fwd.count
for s in samples: print s.pair.rev.count

# Print MD5 checksums #
for s in samples: print s.pair.fwd.md5
for s in samples: print s.pair.rev.md5

################################ Status report ################################
# How far did we run things #
for s in samples: print "Raw:",                s, bool(s.pair)
for s in samples: print "First QC:",           s, bool(s.pair.fwd.fastqc.results)
for s in samples: print "Cleaned:",            s, bool(s.quality_checker.results)
for s in samples: print "Second QC:",          s, bool(s.clean.fwd.fastqc.results)
for s in samples: print "Initial taxa:",       s, bool(s.kraken.results)
for s in samples: print "Solo-assembly:",      s, bool(s.assembly.results)
for s in samples: print "Mono-mapping:",       s, bool(s.mono_mapper.results)
for k,v in proj.assemblies.items(): print "Co-assembly %i:"%k, proj, bool(v.results)
for s in samples: print "Map to co-assembly:", s, bool(s.mapper.results)

################################ Preprocessing ################################
# Clean #
for s in samples:
    print "Cleaning sample '%s'" % s.name
    s.quality_checker.run()

########################## Link from Taito to Sisu ############################
old = "/homeappl/home/bob/"
new = "/wrk/alice/"
for s in samples:
    s.clean.fwd.link_from(s.clean.fwd.path.replace(old, new))
    s.clean.rev.link_from(s.clean.rev.path.replace(old, new))
    s.singletons.link_from(s.singletons.path.replace(old, new))

############################### Co-Assemblies #################################
params = dict(machines=12, cores=12*24, time='12:00:00', partition='small')
proj.runner.run_slurm(steps=['assembly_51.run'], job_name=proj.name+'_ray_51', **params)
proj.runner.run_slurm(steps=['assembly_61.run'], job_name=proj.name+'_ray_61', **params)
proj.runner.run_slurm(steps=['assembly_71.run'], job_name=proj.name+'_ray_71', **params)
proj.runner.run_slurm(steps=['assembly_81.run'], job_name=proj.name+'_ray_81', **params)

################################# Solo-Assembly ###############################
params = dict(steps=['assembly.run'], machines=6, cores=6*24, time='12:00:00', partition='small')
for s in samples: s.runner.run_slurm(job_name = s.name+'_ray', **params)

########################## Link from Sisu to Taito ############################
old = "/homeappl/home/alice/"
new = "/wrk/bob/"
proj.p.assembly_dir.link_from(proj.p.assembly_dir.path.replace(old, new), safe=True)
for s in samples:  s.p.assembly_dir.link_from(s.p.assembly_dir.path.replace(old, new), safe=True)

################################ Merge-Assembly ###############################
params = dict(machines=1, cores=24, time='14-00:00:00', partition='longrun', constraint='hsw', memory=120000)
proj.runner.run_slurm(steps=['merged.run'], job_name="ice_newbler", **params)

################################## Mappings ###################################
params = dict(machines=1, cores=1, threads=24, time='14-00:00:00', partition='longrun', constraint='hsw', memory=120000)
for s in samples: s.runner.run_slurm(steps=['mapper_51.run'],     job_name=s.name + "_co_51_map", **params)
for s in samples: s.runner.run_slurm(steps=['mapper_61.run'],     job_name=s.name + "_co_61_map", **params)
for s in samples: s.runner.run_slurm(steps=['mapper_71.run'],     job_name=s.name + "_co_71_map", **params)
for s in samples: s.runner.run_slurm(steps=['mapper_81.run'],     job_name=s.name + "_co_81_map", **params)
for s in samples: s.runner.run_slurm(steps=['mapper_merged.run'], job_name=s.name + "_merge_map", **params)
for s in samples: s.runner.run_slurm(steps=['mono_mapper.run'],   job_name=s.name + "_merge_map", **params)