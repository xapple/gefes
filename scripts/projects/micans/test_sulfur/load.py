#!/usr/bin/env python2
# -*- coding: utf-8 -*-

"""
A script to load the basic objects for the granular sludge project.
"""

# Built-in modules #
import shutil

# First party modules #
from plumbing.processes import prll_map
from plumbing.timer     import Timer
from plumbing.autopaths import FilePath

# Internal modules #
import gefes

# Third party modules #
import sh

# Home #
import os
home = os.environ.get('HOME', '~') + '/'

###############################################################################
# Parameters #
gefes.groups.samples.Sample.default_cleaner         = "dummy"
gefes.groups.samples.Sample.default_assembler       = "dummy"
gefes.groups.aggregates.Aggregate.default_assembler = "dummy"
gefes.groups.samples.Sample.default_mapper          = "bbmap"

#################################### Load #####################################
# Load project #
proj = gefes.load("~/deploy/gefes/metadata/json/projects/micans/test_sulfur/", raw_files_must_exist=False)