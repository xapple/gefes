#!/usr/bin/env python

"""
A script to build the small test files.
They should contain only 2000 sequences.

Written by Lucas Sinclair.
Kopimi.

You can use this script from the shell like this:
$ ./gen_test.py
"""

# Built-in modules #
import time, datetime

# Internal modules #
from gefes import projects

# Third party modules #
from shell_command import shell_output
import playdoh

# Timer #
now = time.time()

###############################################################################
print "Making test files"

# Do it #
pairs = []
pairs += [(projects['humic'][i].fwd_path, projects['test'][i].fwd_path) for i in range(3)]
pairs += [(projects['humic'][i].rev_path, projects['test'][i].rev_path) for i in range(3)]
process = lambda x : shell_output('zcat %s |head -n 4000| gzip > %s' % (x[0],x[1]))

# Run it in parallel #
playdoh.map(process, pairs, cpu=len(pairs))

# Report Success #
run_time = datetime.timedelta(seconds=round(time.time()-now))
print "\033[0;32mRun time: '%s'\033[0m" % (run_time)