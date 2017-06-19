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
#proj1.merged.results.hit_profile.run()

#################################### Report #####################################
#proj.first.report.generate()
#path = "Dropbox/Granular sludge delivery/projects/granular_sludge_1/reports/samples/s01.pdf"
#proj.first.report.output_path.copy(gefes.home + path)

#################################### Profile #####################################
#proj1.merged.results.hit_profile.run()

#################################### Bundle #####################################
from gefes.distribute.bundle import Bundle
bundle = Bundle("granular_sludge", proj1.samples + proj2.samples)
with Timer(): bundle.run()
