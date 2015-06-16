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
for p in projects: print "Co-assembly:",       p, bool(p.assembly.results)
for s in samples: print "Map to co-assembly:", s, bool(s.mapper.results)

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

########################## Link from Sisu to Taito ############################
old = "/homeappl/home/bob/"
new = "/wrk/alice/"
for s in samples:
    s.clean.fwd.link_from(s.clean.fwd.path.replace(old, new))
    s.clean.rev.link_from(s.clean.rev.path.replace(old, new))
    s.singletons.link_from(s.singletons.path.replace(old, new))

################################# Co-Assembly #################################
for proj in projects:
    proj.runner.run_slurm(steps     = ['assembly.run'],
                          machines  = 24,
                          cores     = 24*24,
                          time      = '12:00:00',
                          partition = 'small',
                          job_name  = proj.name + '_ray71',
                          email     = False)

################################# Solo-Assembly ###############################
params = {'steps'     : ['assembly.run'],
          'machines'  : 8,
          'cores'     : 8*24,
          'time'      : '12:00:00',
          'partition' : 'small'}

for s in samples: s.runner.run_slurm(job_name = s.name+'_ray', **params)