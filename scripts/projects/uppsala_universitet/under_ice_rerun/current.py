#!/usr/bin/env python2

"""
A script to run small snippets of code on the under ice rerun project.
"""

# Modules #
import gefes, shutil

# Three projects #
bt = gefes.load("~/deploy/gefes/metadata/json/projects/uppsala_universitet/bt/")
lb = gefes.load("~/deploy/gefes/metadata/json/projects/uppsala_universitet/lb/")
kt = gefes.load("~/deploy/gefes/metadata/json/projects/uppsala_universitet/kt/")

# Samples combined #
samples = tuple(bt.samples + lb.samples + kt.samples)
projects = (bt, lb, kt)

###############################################################################
print "Hello world"