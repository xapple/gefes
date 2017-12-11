#!/usr/bin/env python2
# -*- coding: utf-8 -*-

"""
A script to load the basic objects for the granular sludge project.
"""

# Built-in modules #
import os, inspect

# Internal modules #
import gefes

# Third party modules #

###############################################################################
# Get current directory #
file_name = inspect.getframeinfo(inspect.currentframe()).filename
this_dir  = os.path.dirname(os.path.abspath(file_name)) + '/'
execfile(this_dir + "load.py")

# Home #
import os
home = os.environ.get('HOME', '~') + '/'

#################################### Load #####################################
# Set up the cleaned reads #
# /GEFES/raw/projects/micans/deep_sulfur/3726_micans/200_cleaned_reads/ONKKR15_S3_L001_1P.gz
# /GEFES/raw/projects/micans/deep_sulfur/3726_micans/200_cleaned_reads/ONKKR15_S3_L001_1U.gz

def uncompress_clean(s):
    # The paired reads #
    dest   = s.clean.fwd
    source = home + "GEFES/raw/projects/micans/deep_sulfur/3726_micans/200_cleaned_reads/%s_S%i_L001_1P.gz"
    source = source % (s.long_name, s.num)
    print source, dest # source.ungzip_to(dest)
    dest   = s.clean.rev
    source = home + "GEFES/raw/projects/micans/deep_sulfur/3726_micans/200_cleaned_reads/%s_S%i_L001_2P.gz"
    source = source % (s.long_name, s.num)
    print source, dest # source.ungzip_to(dest)
    # The singletons #
    dest   = s.singletons
    source = home + "GEFES/raw/projects/micans/deep_sulfur/3726_micans/200_cleaned_reads/%s_S%i_L001_1U.gz"
    source = source % (s.long_name, s.num)
    print source, dest # source.ungzip_to(dest)
    source = home + "GEFES/raw/projects/micans/deep_sulfur/3726_micans/200_cleaned_reads/%s_S%i_L001_2U.gz"
    source = source % (s.long_name, s.num)
    print source, dest # source.ungzip_to(dest, mode='a')

for s in proj: uncompress_clean(s)

