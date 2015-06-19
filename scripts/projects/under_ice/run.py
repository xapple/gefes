#!/usr/bin/env python2

"""
A script to contain the procedure for running the soda evaluation project.
"""

# Built-in modules #

# Internal modules #
import gefes

#################################### Load #####################################
# Three projects #
bt = gefes.projects['under_ice_bt'].load()
lb = gefes.projects['under_ice_lb'].load()
kt = gefes.projects['under_ice_kt'].load()
for s in bt.samples: s.load()
for s in lb.samples: s.load()
for s in kt.samples: s.load()
samples = tuple(bt.samples + lb.samples + kt.samples)
projects = (bt, lb, kt)

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
for s in samples:  print "Raw:",                s, bool(s.pair)
for s in samples:  print "First QC:",           s, bool(s.pair.fwd.fastqc.results)
for s in samples:  print "Cleaned:",            s, bool(s.quality_checker.results)
for s in samples:  print "Second QC:",          s, bool(s.clean.fwd.fastqc.results)
for s in samples:  print "Initial taxa:",       s, bool(s.kraken.results)
for s in samples:  print "Solo-assembly:",      s, bool(s.assembly.results)
for s in samples:  print "Mono-mapping:",       s, bool(s.mono_mapper.results)
for p in projects: print '\n'.join(["Co-assembly %i: %s %s" % (k, p, bool(v.results)) for k,v in p.assemblies.items()])
for s in samples:  print "Map to co-assembly:", s, bool(s.mapper.results)

################################ Preprocessing ################################
# Clean #
for s in bt:
    print "Cleaning sample '%s'" % s.name
    s.quality_checker.run()
for s in lb:
    print "Cleaning sample '%s'" % s.name
    s.quality_checker.run()
for s in kt:
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
params = dict(machines=24, cores=24*24, time='12:00:00', partition='small')
for proj in projects: proj.runner.run_slurm(steps=['assembly_51.run'], job_name=proj.name+'_ray_51', **params)
for proj in projects: proj.runner.run_slurm(steps=['assembly_61.run'], job_name=proj.name+'_ray_61', **params)
for proj in projects: proj.runner.run_slurm(steps=['assembly_71.run'], job_name=proj.name+'_ray_71', **params)
for proj in projects: proj.runner.run_slurm(steps=['assembly_81.run'], job_name=proj.name+'_ray_81', **params)

############################### Solo-Assemblies ###############################
params = dict(steps=['assembly.run'], machines=8, cores=8*24, time='12:00:00', partition='small')
for s in samples: s.runner.run_slurm(job_name = s.name+'_ray', **params)

########################## Link from Sisu to Taito ############################
old = "/homeappl/home/alice/"
new = "/wrk/bob/"
for p in projects:
    print "rsync -av --progress %s %s" % (p.p.assembly_dir.path.replace(old, new), p.p.assembly_dir)
for s in samples:
    print "rsync -av --progress %s %s" % (s.p.assembly_dir.path.replace(old, new), s.p.assembly_dir)

################################ Merge-Assembly ###############################
params = dict(machines=1, cores=1, time='14-00:00:00', partition='longrun',
              memory=120000, constraint='hsw')
proj.runner.run_slurm(steps=['merged.run'], job_name="ice_newbler", **params)

########################### Mappings to Co-Assembly ###########################
params = dict(machines=1, cores=1, time='00:30:00', partition='test',
              threads=6, mem_per_cpu=5300, constraint='hsw')
s = samples[0]
s.runner.run_slurm(steps=[{'mapper_71.run':{'cpus':6}}], job_name=s.name + "_co_71_map_test", **params)

params = dict(machines=1, cores=1, time='14-00:00:00', partition='longrun',
              threads=6, mem_per_cpu=5300, constraint='hsw')
for s in samples: s.runner.run_slurm(steps=['mapper_51.run'],     job_name=s.name + "_co_51_map", **params)
for s in samples: s.runner.run_slurm(steps=['mapper_61.run'],     job_name=s.name + "_co_61_map", **params)
for s in samples[1:]: s.runner.run_slurm(steps=['mapper_71.run'],     job_name=s.name + "_co_71_map", **params)
for s in samples: s.runner.run_slurm(steps=['mapper_81.run'],     job_name=s.name + "_co_81_map", **params)
for s in samples: s.runner.run_slurm(steps=['mapper_merged.run'], job_name=s.name + "_merge_map", **params)
for s in samples: s.runner.run_slurm(steps=['mono_mapper.run'],   job_name=s.name + "_merge_map", **params)