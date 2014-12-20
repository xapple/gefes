#!/usr/bin/env python2

"""
A script to build three small test files for
running the pipeline quickly when developing.
They should contain about 2000 sequences each.

Written by Lucas Sinclair.
Kopimi.

You can use this script from the shell like this:
$ ./generate_test.py
"""

# Built-in modules #
import time, datetime

# Internal modules #
from gefes import projects

# First party modules #
from plumbing.color import Color

###############################################################################
# Message #
now = time.time()
print "Making test files"

# Choose a project #
proj = projects['soda_eval']
proj.load()

# Do it #
projects['test'].load()
source_to_dest = []
source_to_dest += [(proj[i].pair.fwd, projects['test'][i].pair.fwd) for i in (0,1,2)]
source_to_dest += [(proj[i].pair.rev, projects['test'][i].pair.rev) for i in (0,1,2)]

# Downsample #
def downsample(source, dest): dest.write(source[0:2000])

# Run each one #
for source, dest in source_to_dest: downsample(source, dest)

# Report Success #
run_time = datetime.timedelta(seconds=round(time.time()-now))
print Color.grn + ("Run time: '%s'" % run_time) + Color.end

# Count sequences #
for source, dest in source_to_dest: print dest.prefix + ': ' + str(len(dest)) + ' sequences'