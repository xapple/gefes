#!/usr/bin/env python2
# -*- coding: utf-8 -*-

"""
A script to load the basic objects for the granular sludge project.
"""

# Built-in modules #

# First party modules #
from plumbing.processes import prll_map
from plumbing.timer import Timer

# Internal modules #
import gefes

#################################### Load #####################################
# Load two projects #
proj1 = gefes.load("~/deploy/gefes/metadata/json/projects/epfl/granular_sludge_1/")
proj2 = gefes.load("~/deploy/gefes/metadata/json/projects/epfl/granular_sludge_2/")
projects = (proj1, proj2)

# Parameters #
for s in proj1: s.default_cleaner = 'window'
for s in proj2: s.default_cleaner = 'window'
