#!/usr/bin/env python2
# -*- coding: utf-8 -*-

"""
A script to run small snippets of code on the granular sludge project.
"""

# Built-in modules #
import os, inspect

# First party modules #
from plumbing.processes import prll_map
from plumbing.timer import Timer

# Internal modules #
import gefes

#################################### Load #####################################
# Get current directory #
file_name = inspect.getframeinfo(inspect.currentframe()).filename
this_dir  = os.path.dirname(os.path.abspath(file_name)) + '/'
execfile(this_dir + "load.py")

# Home #
import os
home = os.environ.get('HOME', '~') + '/'

###############################################################################
source = home + "/GEFES/raw/projects/micans/deep_sulfur/3726_micans/300_assemblies/320_coassemblies/322_metaspades/contigs.fasta"
proj.merged.p.contigs.link_from(source)
