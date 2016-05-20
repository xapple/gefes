#!/usr/bin/env python2

"""
A script to load the main objects for the soda rerun project.
"""

# Built-in modules #

# Internal modules #
import gefes

#################################### Load #####################################
# One project #
proj = gefes.projects['soda_rerun'].load()
samples = proj.samples
for s in samples: s.load()
