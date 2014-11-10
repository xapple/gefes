#!/usr/bin/env python2

"""
A script to setup some symbolic links on the Quince server for the strain mock project.
"""

# Modules #
import os, gefes
from plumbing.autopaths import DirectoryPath

###############################################################################
home = os.environ['HOME'] + '/'
result_dir = DirectoryPath(home + "storage/ino/projects/2014hack/StrainMock/out_41_cutup_10K_nodup/map/")
p = gefes.projects['strain_mock'].load()
for sample in p:
    sample.load()
    computed_dir = result_dir + sample.name + ".fasta.gz/bowtie2/"
    pairs = [
        (computed_dir + "asm_pair-smd.metrics",   sample.mapper.p.metrics),
        (computed_dir + "asm_pair-smds.bam",      sample.mapper.p.map_smds_bam),
        (computed_dir + "asm_pair-smds.bam.bai",  sample.mapper.p.map_smds_bai),
        (computed_dir + "asm_pair-smds.coverage", sample.mapper.p.coverage),
    ]
    for orig,dest in pairs: os.symlink(orig, dest) # print orig+'\n'+dest+'\n\n'