#!/usr/bin/env python2

"""
A script to build three small test files.
They should contain only 2000 sequences each.

Written by Lucas Sinclair.
Kopimi.

You can use this script from the shell like this:
$ ./generate_test.py
"""

# Built-in modules #
import time, datetime

# Internal modules #
from plumbing.color import Color
from gefes import projects

# Third party modules #
from shell_command import shell_output
import playdoh

###############################################################################
# Message #
now = time.time()
print "Making test files"

# Do it #
source_to_dest = []
source_to_dest += [(projects['alinen'][i].pair.fwd, projects['test'][i].pair.fwd) for i in (0,1,2)]
source_to_dest += [(projects['alinen'][i].pair.rev, projects['test'][i].pair.rev) for i in (0,1,2)]
downsample = lambda x : shell_output('zcat %s |head -n 4000| gzip > %s' % (x[0],x[1]))

# Run it in parallel #
playdoh.map(downsample, source_to_dest, cpu=len(downsample))

# Report Success #
run_time = datetime.timedelta(seconds=round(time.time()-now))
print Color.grn + ("Run time: '%s'" % run_time) + Color.end

# Count sequences #
for source, dest in source_to_dest: print dest.prefix + ': ' + len(dest) + ' sequences'