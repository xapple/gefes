#!/usr/bin/env python2
# -*- coding: utf-8 -*-

"""
A script to load the basic objects for the granular sludge project.
"""

# Built-in modules #
import os, inspect

# First party modules #
from plumbing.autopaths import FilePath

# Internal modules #
import gefes

# Third party modules #

#################################### Load #####################################
# Get current directory #
file_name = inspect.getframeinfo(inspect.currentframe()).filename
this_dir  = os.path.dirname(os.path.abspath(file_name)) + '/'
execfile(this_dir + "load.py")

# Home #
import os
home = os.environ.get('HOME', '~') + '/'

#################################### Generate #################################
# Real project #
real_proj = gefes.load("~/deploy/gefes/metadata/json/projects/micans/deep_sulfur/")

# Do it #
source_to_dest  = []
source_to_dest += [(real_proj[i].clean.fwd, proj[i].clean.fwd) for i in (1,2,3)]
source_to_dest += [(real_proj[i].clean.rev, proj[i].clean.rev) for i in (1,2,3)]

# Downsample #
def downsample(source, dest): dest.write(itertools.islice(source, 0, 2000))

# Run each one #
for source, dest in source_to_dest: downsample(source, dest)

# Count sequences #
for source, dest in source_to_dest: print dest.prefix + ': ' + str(len(dest)) + ' sequences'
