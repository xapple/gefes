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
import time, datetime, itertools

# Internal modules #
from gefes import projects

# First party modules #
from plumbing.color import Color
from fasta import FASTA

###############################################################################
# Message #
now = time.time()
print "Making test files"

# Choose a project #
proj = projects['alinen']
proj.load()

# Load stuff #
projects['test'].load()
for s in projects['test']: s.load()

# Do it #
source_to_dest = []
source_to_dest += [(proj[i].pair.fwd, projects['test'][i].pair.fwd) for i in (0,1,2)]
source_to_dest += [(proj[i].pair.rev, projects['test'][i].pair.rev) for i in (0,1,2)]

# Downsample #
def downsample(source, dest): dest.write(itertools.islice(source, 0, 2000))

# Run each one #
for source, dest in source_to_dest: downsample(source, dest)

# Count sequences #
for source, dest in source_to_dest: print dest.prefix + ': ' + str(len(dest)) + ' sequences'

# We need some contigs too... #
contigs = FASTA(proj.base_dir + "assembly/ray/71/output/Contigs.fasta").parse()
destinations = [FASTA(projects['test'][i].assembly.p.Contigs) for i in (0,1,2)]
for dest in destinations: dest.write([contigs.next() for x in range(10)])

# Count sequences #
for dest in destinations: print dest.prefix + ': ' + str(len(dest)) + ' sequences'

# Report Success #
run_time = datetime.timedelta(seconds=round(time.time()-now))
print Color.grn + ("Run time: '%s'" % run_time) + Color.end