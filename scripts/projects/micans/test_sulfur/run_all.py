#!/usr/bin/env python2
# -*- coding: utf-8 -*-

"""
A script to run all the procedure on the granular sludge project.

First create the JSON files from the Excel file.
"""

# Built-in modules #
import os, shutil

# Internal modules #
import gefes

# First party modules #
from plumbing.processes import prll_map
from plumbing.timer     import Timer
from plumbing.autopaths import FilePath

# Third party modules #
import sh
from tqdm import tqdm

#################################### Load #####################################
# Parameters #
gefes.groups.samples.Sample.default_cleaner         = "dummy"
gefes.groups.samples.Sample.default_assembler       = "dummy"
gefes.groups.aggregates.Aggregate.default_assembler = "dummy"
gefes.groups.samples.Sample.default_mapper          = "bbmap"

# Load project #
proj = gefes.load("~/deploy/gefes/metadata/json/projects/micans/test_sulfur/",
                  raw_files_must_exist=False)


###############################################################################
################################ Print status #################################
###############################################################################
proj.status.print_long()
proj.status.print_short()

################################## Validate ###################################
with Timer(): prll_map(lambda s: s.clean.validate(), proj.samples, 5)    # xhxx

################################## FASTQC #####################################
for s in tqdm(proj):                                                     # xhxx
    print "\n FastQC on sample '%s'" % s.name
    s.pair.fwd.fastqc.run(cpus=4)
    s.pair.rev.fastqc.run(cpus=4)
    s.clean.fwd.fastqc.run(cpus=4)
    s.clean.rev.fastqc.run(cpus=4)

################################ Co-mappings ##################################
print("# Co-mapping #")                                                  # xhxx
with Timer(): proj[1].mapper_merged.run(cpus=32)
with Timer(): proj[2].mapper_merged.run(cpus=32)
with Timer(): proj[3].mapper_merged.run(cpus=32)

################################ Mono-mappings ################################
proj[1].mono_mapper.run(cpus=16)                                         # xhxx
proj[2].mono_mapper.run(cpus=16)
proj[3].mono_mapper.run(cpus=16)

################################# Binning #####################################
with Timer(): proj.merged.results.binner.run()                           # xhxx


###############################################################################
################################## Analysis ###################################
###############################################################################

################################## CheckM #####################################
for b in tqdm(proj.merged.results.binner.results.bins): b.evaluation.run(cpus=32)

################################ Prodigal #####################################
for c in tqdm(proj.merged.results.contigs): c.proteins.run()
for b in tqdm(proj.merged.results.binner.results.bins): b.faa

################################ Phylophlan ###################################
proj.merged.results.binner.results.taxonomy.run(cpus=32)

################################ Pfam #########################################
bins = proj.merged.results.binner.results.bins
with Timer(): prll_map(lambda b: b.pfams.run(cpus=1), bins, 45)          # xhxx

################################ Profile ######################################
with Timer(): proj.merged.results.hit_profile.run()

################################ Phylosift ####################################
#for c in proj.merged.results.contigs: c.taxonomy.run()


###############################################################################
################################## Reports ####################################
###############################################################################

################################## Project ####################################
proj.report.generate()

################################## Samples ####################################
for s in tqdm(proj): # xhxx
    print "Report on sample '%s'" % s.name
    s.report.generate()
with Timer(): prll_map(lambda s: s.report.generate(), proj, 32)          # xhxx

################################# Assembly ####################################
proj.merged.results.report.generate()

#################################### Bins #####################################
bins = proj.merged.results.binner.results.bins
with Timer(): prll_map(lambda b: b.report.generate(), bins, 32)


###############################################################################
################################## Delivery ###################################
###############################################################################

################################## Bundle #####################################
from gefes.distribute.bundle import Bundle
bundle = Bundle("deep_sulfur", proj.samples)
with Timer(): bundle.run()

################################ Extra files ##################################
path = FilePath(gefes.home + "deploy/gefes/metadata/excel/projects/micans/deep_sulfur/metadata.xlsx")
shutil.copy(path, bundle.p.samples_xlsx)

################################## Upload #####################################
from gefes.distribute.dropbox import DropBoxRclone
dbx_sync = DropBoxRclone(bundle.base_dir, '/Deep sulfur delivery')
with Timer(): dbx_sync.run()
print("Total delivery: %s" % bundle.base_dir.size)
