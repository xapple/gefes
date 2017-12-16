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

########################### Set up the cleaned reads ##########################
# /GEFES/raw/projects/micans/deep_sulfur/3726_micans/200_cleaned_reads/ONKKR15_S3_L001_1P.gz
# /GEFES/raw/projects/micans/deep_sulfur/3726_micans/200_cleaned_reads/ONKKR15_S3_L001_1U.gz

def uncompress_clean(s):
    # The paired reads #
    dest   = s.clean.fwd
    source = home + "GEFES/raw/projects/micans/deep_sulfur/3726_micans/200_cleaned_reads/%s_S%i_L001_1P.gz"
    source = FilePath(source % (s.long_name, s.num))
    source.ungzip_to(dest)
    dest   = s.clean.rev
    source = home + "GEFES/raw/projects/micans/deep_sulfur/3726_micans/200_cleaned_reads/%s_S%i_L001_2P.gz"
    source = FilePath(source % (s.long_name, s.num))
    source.ungzip_to(dest)
    # The singletons #
    dest   = s.singletons
    source = home + "GEFES/raw/projects/micans/deep_sulfur/3726_micans/200_cleaned_reads/%s_S%i_L001_1U.gz"
    source = FilePath(source % (s.long_name, s.num))
    source.ungzip_to(dest)
    source = home + "GEFES/raw/projects/micans/deep_sulfur/3726_micans/200_cleaned_reads/%s_S%i_L001_2U.gz"
    source = FilePath(source % (s.long_name, s.num))
    source.ungzip_to(dest, mode='a')

for s in proj: uncompress_clean(s)

######################### Set up the mono-assemblies ###########################
# /GEFES/raw/projects/micans/deep_sulfur/3726_micans/300_assemblies/310_single_samples/313_metaspades/ONKPVA08_S2/contigs.fasta

def setup_monoassembly(s):
    source = home + "GEFES/raw/projects/micans/deep_sulfur/3726_micans/300_assemblies/310_single_samples/313_metaspades/%s_S%i/contigs.fasta"
    source = FilePath(source % (s.long_name, s.num))
    s.assembly.p.contigs.link_from(source)
    s.assembly.run()

for s in proj: setup_monoassembly(s)

########################### Set up the co-assembly #############################
source = home + "GEFES/raw/projects/micans/deep_sulfur/3726_micans/300_assemblies/320_coassemblies/322_metaspades/contigs.fasta"
proj.merged.p.contigs.link_from(source)
source = home + "GEFES/raw/projects/micans/deep_sulfur/3726_micans/300_assemblies/301_all_assemblies/spades_coassembly.le2500bp.contigs.fa"
proj.merged.p.filtered.link_from(source)
