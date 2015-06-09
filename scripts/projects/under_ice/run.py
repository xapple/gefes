#!/usr/bin/env python2

"""
A script to contain the procedure for running the soda evaluation project.
"""

# Built-in modules #

# Internal modules #
import gefes

# Constants #
bt = gefes.projects['under_ice_bt'].load()
lb = gefes.projects['under_ice_lb'].load()
kt = gefes.projects['under_ice_kt'].load()
for s in bt.samples: s.load()
for s in lb.samples: s.load()
for s in kt.samples: s.load()
samples = bt.samples + lb.samples + kt.samples

################################ Preprocessing ################################
# Print number of sequences #
for s in samples: print s.pair.fwd.count
for s in samples: print s.pair.rev.count

# Print MD5 checksums #
for s in samples: print s.pair.fwd.md5
for s in samples: print s.pair.rev.md5

# Clean #
for s in samples: s.runner.run_slurm()