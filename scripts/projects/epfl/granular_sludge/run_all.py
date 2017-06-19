#!/usr/bin/env python2
# -*- coding: utf-8 -*-

"""
A script to run all the procedure on the granular sludge project.
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
# Load two projects #
proj1 = gefes.load("~/deploy/gefes/metadata/json/projects/epfl/granular_sludge_1/")
proj2 = gefes.load("~/deploy/gefes/metadata/json/projects/epfl/granular_sludge_2/")
projects = (proj1, proj2)
proj = proj1 + proj2

###############################################################################
print("# Get information for excel file #")
for s in proj: print s.short_name
for s in proj: print s.pair.fwd.count
for s in proj: print s.pair.rev.count
for s in proj: print s.pair.fwd.md5
for s in proj: print s.pair.rev.md5
for s in proj: print len(s.pair.fwd.first)

print("# Get special sample name for proj2 #")
from collections import OrderedDict
d = open("/home/lucas/GEFES//raw/projects/epfl/granular_sludge/2/correspondances.csv").read()
d = OrderedDict((l.split()[0][:-13] + '.fastq.gz', l.split()[1]) for l in d.split('\n') if l)
for s in proj2: print s.short_name
for s in proj2: print d[s.pair.fwd.filename]

print("# Convert Old 1.3 PHRED format for proj2 #")
with Timer(): prll_map(lambda s: s.pair.fwd.phred_13_to_18(), proj2)
with Timer(): prll_map(lambda s: s.pair.rev.phred_13_to_18(), proj2)

print("# Check PHRED formats #")
for s in proj1: assert s.pair.fwd.guess_phred_format() == 'Sanger'
for s in proj2: assert s.pair.rev.guess_phred_format() == 'Sanger'


###############################################################################
################################ Print status #################################
###############################################################################
proj1.status.print_long()
proj2.status.print_short()

################################ Preprocessing ################################
print("# Starting cleaning of samples #")
with Timer(): prll_map(lambda s: s.quality_checker.run(), proj1)
with Timer(): prll_map(lambda s: s.quality_checker.run(), proj2)

################################ Co-assembly ##################################
print("# Co-assembly #")
with Timer(): proj1.merged.run() # 40 minutes
with Timer(): proj2.merged.run() # 1 hour 26 minutes

################################ Co-mappings ###################################
print("# Co-mapping #")
with Timer(): prll_map(lambda s: s.mapper.run(cpus=4), proj1, cpus=16) # x minutes
with Timer(): prll_map(lambda s: s.mapper.run(cpus=4), proj2, cpus=16) # x minutes
for s in tqdm(proj): s.mapper.run(cpus=16) # 3 hours

################################# Binning #####################################
with Timer(): proj1.merged.results.binner.run() # 3.5 hours
with Timer(): proj2.merged.results.binner.run() # 11.5 hours

################################## FASTQC #####################################
for s in tqdm(proj): # 1.5 hours
    print "\n FastQC on sample '%s'" % s.name
    s.pair.fwd.fastqc.run(cpus=4)
    s.pair.rev.fastqc.run(cpus=4)
    s.clean.fwd.fastqc.run(cpus=4)
    s.clean.rev.fastqc.run(cpus=4)

################################## Kraken #####################################
for s in tqdm(proj): # Skipped
    print "\n Kraken on sample '%s'" % s.name
    s.kraken.run(cpus=4)

################################ Mono-assemblies ##############################
for s in tqdm(proj): # xx hours   (NOTE: sample s27 empty)
    print "\n Mono assembly on sample '%s'" % s.name
    s.assembly.run(cpus=32)

################################ Mono-mappings ################################
for s in tqdm(proj): # xx hours
    print "\n Mono mapping on sample '%s'" % s.name
    s.mono_mapper.run(cpus=32)


###############################################################################
################################## Analysis ###################################
###############################################################################

################################## CheckM #####################################
for b in tqdm(proj1.merged.results.binner.results.bins): b.evaluation.run(cpus=32)
for b in tqdm(proj2.merged.results.binner.results.bins): b.evaluation.run(cpus=32)

################################ Prodigal #####################################
for c in tqdm(proj1.merged.results.contigs): c.proteins.run()
for c in tqdm(proj2.merged.results.contigs): c.proteins.run()
for b in tqdm(proj1.merged.results.binner.results.bins): b.faa
for b in tqdm(proj2.merged.results.binner.results.bins): b.faa

################################ Phylophlan ###################################
proj1.merged.results.binner.results.taxonomy.run(cpus=32)
proj2.merged.results.binner.results.taxonomy.run(cpus=32)

################################ Pfam #########################################
bins = proj1.merged.results.binner.results.bins
with Timer(): prll_map(lambda b: b.pfams.run(cpus=1), bins, 45)          # 0h12
bins = proj2.merged.results.binner.results.bins
with Timer(): prll_map(lambda b: b.pfams.run(cpus=1), bins, 45)          # 0h22

################################ Profile ######################################
with Timer(): proj1.merged.results.hit_profile.run()
with Timer(): proj2.merged.results.hit_profile.run()

################################ Phylosift ####################################
#contigs = proj1.merged.results.contigs + proj2.merged.results.contigs
#for c in proj1.merged.results.contigs: c.taxonomy.run()
#for c in proj2.merged.results.contigs: c.taxonomy.run()


###############################################################################
################################## Reports ####################################
###############################################################################

################################## Projects ###################################
proj1.report.generate()
proj2.report.generate()

################################## Samples ####################################
for s in tqdm(proj): # 1h43
    print "Report on sample '%s'" % s.name
    s.report.generate()
with Timer(): prll_map(lambda s: s.report.generate(), proj, 32) # 1h43

################################# Assemblies ##################################
proj1.merged.results.report.generate()
proj2.merged.results.report.generate()

#################################### Bins #####################################
bins = proj1.merged.results.binner.results.bins + proj2.merged.results.binner.results.bins
with Timer(): prll_map(lambda b: b.report.generate(), bins, 32)


###############################################################################
################################## Delivery ###################################
###############################################################################

################################## Bundle #####################################
from gefes.distribute.bundle import Bundle
bundle = Bundle("granular_sludge", proj1.samples + proj2.samples)
with Timer(): bundle.run()

################################ Extra files ##################################
path = gefes.home + "deploy/gefes/metadata/excel/projects/epfl/granular_sludge/metadata.xlsx"
shutil.copy(path, bundle.p.samples_xlsx)

################################## Upload #####################################
from gefes.distribute.dropbox import DropBoxRclone
dbx_sync = DropBoxRclone(bundle.base_dir, '/Granular sludge delivery')
with Timer(): dbx_sync.run()
print("Total delivery: %s" % bundle.base_dir.size)
