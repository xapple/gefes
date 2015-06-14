#!/usr/bin/env python2

"""
A script to contain the procedure for running the soda evaluation project.
"""

# Built-in modules #

# Internal modules #
import gefes

# Third party modules #
from tqdm import tqdm

#################################### Load #####################################
# One project #
proj = gefes.projects['soda_rerun'].load()
samples = proj.samples
for s in samples: s.load()

############################## Special Merging ################################
for s in tqdm(samples): s.merge_lanes()

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
print                   "Co-assembly:",     proj, bool(s.assembly.results)
for s in samples: print "Map to co-assembly:", s, bool(s.mapper.results)

################################ Preprocessing ################################
# Clean #
for s in samples:
    print "Cleaning sample '%s'" % s.name
    s.quality_checker.run()

########################## Link from Sisu to Taito ############################
old = "/homeappl/home/bob/"
new = "/wrk/alice/"
for s in samples:
    print s.clean.fwd.link_from(s.clean.fwd.path.replace(old, new))
    print s.clean.rev.link_from(s.clean.rev.path.replace(old, new))

################################# Co-Assembly #################################
proj.runner.run_slurm(steps     = ['assembly.run'],
                      machines  = 42,
                      cores     = 42*24,
                      time      = '36:00:00',
                      partition = 'large',
                      job_name  = 'alinen_ray',
                      email     = False)
