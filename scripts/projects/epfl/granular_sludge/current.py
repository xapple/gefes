#!/usr/bin/env python2
# -*- coding: utf-8 -*-

"""
A script to run small snippets of code on the granular sludge project.
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
proj = proj1 + proj2

###############################################################################
################################### Stuff #####################################
###############################################################################
#print proj1.merged.report.generate()
#print proj.first.report.generate()
proj1.merged.results.hit_profile.run()