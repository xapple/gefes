#!/usr/bin/env python2

"""
A script to contain the procedure for running the under_ice rerun project.

ipython -i ~/repos/gefes/scripts/projects/under_ice_rerun/run.py
"""

# Built-in modules #
import os

# Internal modules #
import gefes

# Third party modules #
from tqdm import tqdm

# Constants #
user = os.environ.get('USER')

#################################### Load #####################################
# Three projects #
bt = gefes.projects['under_ice_rerun_bt'].load()
lb = gefes.projects['under_ice_rerun_lb'].load()
kt = gefes.projects['under_ice_rerun_kt'].load()
for s in bt.samples: s.load()
for s in lb.samples: s.load()
for s in kt.samples: s.load()
samples = tuple(bt.samples + lb.samples + kt.samples)
projects = (bt, lb, kt)

# Don't run it #
import sys
sys.exit("Copy paste the commands you want in ipython, don't run this script.")

############################## Special Merging ################################
for s in tqdm(bt.samples): s.merge_lanes(remove_orig=True)
for s in tqdm(lb.samples): s.merge_lanes(remove_orig=True)
for s in tqdm(kt.samples): s.merge_lanes(remove_orig=True)

################################## Meta-data ##################################
# Special merging #
# Print number of sequences #
for s in samples:
    s.pair.fwd.__cache__.pop('count')
    print s.pair.fwd.count

# Print MD5 checksums #
for s in samples: print s.pair.fwd.md5
for s in samples: print s.pair.rev.md5

################################ Status report ################################
# How far did we run things #
for s in samples:  print "Raw:",                s, bool(s.pair)
for s in samples:  print "First QC:",           s, bool(s.pair.fwd.fastqc)
for s in samples:  print "Cleaned:",            s, bool(s.quality_checker.results)
for s in samples:  print "Second QC:",          s, bool(s.clean.fwd.fastqc)
for s in samples:  print "Initial taxa:",       s, bool(s.kraken)
for s in samples:  print "Mono-assembly:",      s, bool(s.assembly)
for p in projects: print '\n'.join(["Co-assembly %i: %s %s" % (k, p, bool(v)) for k,v in p.assemblies.items()])
for s in samples:  print "Mono-mapping:",       s, bool(s.mono_mapper.results)
for p in projects: print "Merged assembly:",    p, bool(p.merged.results)
for s,a,m in ((s,a,m) for s in samples for a,m in s.mappers.items()): print "Map %s to %s:"%(s,a), bool(m)
for p,n,a in ((p,n,a) for p in projects for n,a in p.assemblies.items()):
                   print "Binning %s %i:"%(p,n), bool(a.results.binner.p.clustering)
for p in projects: print "Merged binning:",     p, bool(p.merged.results.binner.p.clustering)

################################# Search logs ##################################
from plumbing.common import tail
from plumbing.autopaths import FilePath
for s in samples:
    for log in s.runner.logs:
        out = FilePath(log + 'run.out')
        if out.exists and 'mapper_71.run' in out.contents:
            print '-'*50 + '\n'  + tail(log + 'run.out')

################################ Preprocessing ################################
# Clean #
for s in tqdm(samples):
    print "Starting cleaning of sample '%s'" % s.name
    s.quality_checker.run()

################################## FASTQC #####################################
for s in tqdm(samples):
    print "FastQC on sample '%s'" % s.name
    s.pair.fwd.fastqc.run(cpus=4)
    s.pair.rev.fastqc.run(cpus=4)
    s.clean.fwd.fastqc.run(cpus=4)
    s.clean.rev.fastqc.run(cpus=4)

################################## Kraken #####################################
for s in tqdm(samples):
    print "Kraken on sample '%s'" % s.name
    s.kraken.run(cpus=4)

########################## Link from Taito to Sisu ############################
old = "/homeappl/home/bob/"
new = "/wrk/alice/"
for s in samples:
    s.clean.fwd.link_from(s.clean.fwd.path.replace(old, new))
    s.clean.rev.link_from(s.clean.rev.path.replace(old, new))
    s.singletons.link_from(s.singletons.path.replace(old, new))

############################### Co-Assemblies #################################
params = dict(machines=42, cores=42*24, time='36:00:00', partition='large')
for proj in projects: proj.runner.run_slurm(steps=['assembly_51.run'], job_name=proj.name+'_ray_51', **params)
for proj in projects: proj.runner.run_slurm(steps=['assembly_61.run'], job_name=proj.name+'_ray_61', **params)
for proj in projects: proj.runner.run_slurm(steps=['assembly_71.run'], job_name=proj.name+'_ray_71', **params)
for proj in projects: proj.runner.run_slurm(steps=['assembly_81.run'], job_name=proj.name+'_ray_81', **params)

############################### Solo-Assemblies ###############################
params = dict(steps=['assembly.run'], machines=12, cores=12*24, time='12:00:00', partition='small')
for s in samples: s.runner.run_slurm(job_name = s.name+'_ray', **params)

########################## Link from Sisu to Taito ############################
old = "/homeappl/home/eiler/"
new = "/wrk/lsinclai/"
for p in projects:
    print "rsync -av --progress %s %s" % (p.p.assembly_dir.path.replace(old, new), p.p.assembly_dir)
for s in samples:
    print "rsync -av --progress %s %s" % (s.p.assembly_dir.path.replace(old, new), s.p.assembly_dir)

############################### Merged-Assembly ###############################
bt.merged.run(cpus=4)
lb.merged.run(cpus=4)
kt.merged.run(cpus=4)

################################ All Mappings #################################
params = dict(machines=1, cores=1, time='7-00:00:00', partition='longrun',
              threads=6, mem_per_cpu=5300, constraint='hsw')
for s in samples: s.runner.run_slurm(steps=[{'mapper_51.run':{'cpus':6}}],     job_name=s.name + "_co_51_map", **params)
for s in samples: s.runner.run_slurm(steps=[{'mapper_61.run':{'cpus':6}}],     job_name=s.name + "_co_61_map", **params)
for s in samples: s.runner.run_slurm(steps=[{'mapper_71.run':{'cpus':6}}],     job_name=s.name + "_co_71_map", **params)
for s in samples: s.runner.run_slurm(steps=[{'mapper_81.run':{'cpus':6}}],     job_name=s.name + "_co_81_map", **params)
for s in samples: s.runner.run_slurm(steps=[{'mapper_merged.run':{'cpus':6}}], job_name=s.name + "_merge_map", **params)

params = dict(machines=1, cores=1, time='3-00:00:00', partition='serial',
              threads=6, mem_per_cpu=5300, constraint='hsw')
for s in samples: s.runner.run_slurm(steps=[{'mono_mapper.run':{'cpus':6}}],   job_name=s.name + "_mono_map",  **params)

################################# Binnings ####################################
params = dict(machines=1, cores=1, time='7-00:00:00', partition='longrun',
              threads=12, mem_per_cpu=5300, constraint='hsw')
for p in projects: p.runner.run_slurm(steps=['assembly_51.results.binner.run'], job_name=p.name+'_bin_51', **params)
for p in projects: p.runner.run_slurm(steps=['assembly_61.results.binner.run'], job_name=p.name+'_bin_61', **params)
for p in projects: p.runner.run_slurm(steps=['assembly_71.results.binner.run'], job_name=p.name+'_bin_71', **params)
for p in projects: p.runner.run_slurm(steps=['assembly_81.results.binner.run'], job_name=p.name+'_bin_81', **params)
for p in projects: p.runner.run_slurm(steps=['merged.results.binner.run'],      job_name=p.name+'_bin_04', **params)

################################## Prokka #####################################
for a in p.assemblies.values():
    print "Prokka for project '%s', assembly '%s'" % (p.name, a)
    for c in tqdm(a.results.contigs): c.annotation.run()

################################ Phylosift ####################################
for c in proj.assembly.results.contigs: c.taxonomy.run()

################################## CheckM #####################################
params = dict(machines=1, cores=1, memory=124000, time='1-00:00:00', partition='serial', constraint='hsw')
proj.runner.run_slurm(steps=['merged.results.binner.results.run_all_bin_eval'], job_name="checkm_merged", **params)
for b in tqdm(proj.merged.results.binner.results.bins): b.evaluation.run(cpus=4)

################################## Plots ######################################
for s in tqdm(samples):
    print "Plots for sample '%s'" % s.name
    s.clean.fwd.graphs.length_dist.plot(x_log=False, y_log=True)
    s.clean.rev.graphs.length_dist.plot(x_log=False, y_log=True)
    s.singletons.graphs.length_dist.plot(x_log=False, y_log=True)
    s.assembly.results.contigs_fasta.graphs.length_dist.plot(x_log=True, y_log=True)
for p in tqdm(projects):
    print "Plots for project '%s'" % p.name
    for a in p.assemblies.values(): a.results.contigs_fasta.graphs.length_dist.plot(x_log=True, y_log=True)
    p.merged.results.contigs_fasta.graphs.length_dist.plot(x_log=True, y_log=True)

################################## Report #####################################
for s in tqdm(samples):
    print "Report on sample '%s'" % s.name
    s.report.generate()
for p in projects:
    for a in p.assemblies.values():
        print "Report for project '%s', assembly '%s'" % (proj.name, a)
        a.report.generate()