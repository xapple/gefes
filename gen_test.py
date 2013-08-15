#!/usr/bin/env python

"""
A script to build small test files.

Written by Lucas Sinclair.
Kopimi.

You can use this script from the shell like this:
$ gen.py
"""

# Built-in modules #
import time, datetime

# Internal modules #
from hiseq import projects

# Third party modules #
from shell_command import shell_output

# Timer #
now = time.time()

###############################################################################
# Message #
print "Making test files"
for i in range(3):
    shell_output('zcat %s |head -n 1000000| gzip > %s' % (projects['humic'][i].fwd_path, projects['test'][i].fwd_path))
    shell_output('zcat %s |head -n 1000000| gzip > %s' % (projects['humic'][i].rev_path, projects['test'][i].rev_path))

###############################################################################
run_time = datetime.timedelta(seconds=round(time.time()-now))
print "\033[0;32mRun time: '%s'\033[0m" % (run_time)