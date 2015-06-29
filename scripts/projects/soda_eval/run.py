#!/usr/bin/env python2

"""
A script to contain the procedure for running the soda evaluation project.
"""

# Built-in modules #

# Internal modules #
import gefes

# Third party modules #
from tqdm import tqdm

#################################### Load #####################################
# One project #
proj = gefes.projects['soda_eval'].load()
samples = proj.samples
for s in samples: s.load()

# Don't run it #
import sys
sys.exit("Copy paste the commands you want in ipython, don't run this script.")

################################## Meta-data ##################################
# Print number of sequences #
for s in samples: print s.pair.fwd.count
for s in samples: print s.pair.rev.count

# Print MD5 checksums #
for s in samples: print s.pair.fwd.md5
for s in samples: print s.pair.rev.md5

################################ Status report ################################
# How far did we run things #
for s in samples: print "Raw:",                s, bool(s.pair)
for s in samples: print "First QC:",           s, bool(s.pair.fwd.fastqc.results)
for s in samples: print "Cleaned:",            s, bool(s.quality_checker.results)
for s in samples: print "Second QC:",          s, bool(s.clean.fwd.fastqc.results)
for s in samples: print "Initial taxa:",       s, bool(s.kraken.results)
for s in samples: print "Mono-assembly:",      s, bool(s.assembly.results)
for k,v in proj.assemblies.items(): print "Co-assembly %i:"%k, proj, bool(v.results)
for s in samples: print "Mono-mapping:",       s, bool(s.mono_mapper.results)
print                   "Merged assembly:", proj, bool(proj.merged.results)
for s,a,m in ((s,a,m) for a,m in s.mappers.items() for s in samples): print "Map %s to %s:"%(s,a), bool(m.p.coverage)
for k,v in proj.assemblies.items(): print "Binning %i:"%k, proj, bool(v.results.binner.results)
print                   "Merged binning:", bool(proj.merged.results.binner.results)

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
for s in samples:
    print "Cleaning sample '%s'" % s.name
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
params = dict(machines=12, cores=12*24, time='12:00:00', partition='small')
proj.runner.run_slurm(steps=['assembly_51.run'], job_name=proj.name+'_ray_51', **params)
proj.runner.run_slurm(steps=['assembly_61.run'], job_name=proj.name+'_ray_61', **params)
proj.runner.run_slurm(steps=['assembly_71.run'], job_name=proj.name+'_ray_71', **params)
proj.runner.run_slurm(steps=['assembly_81.run'], job_name=proj.name+'_ray_81', **params)

################################# Solo-Assembly ###############################
params = dict(steps=['assembly.run'], machines=6, cores=6*24, time='12:00:00', partition='small')
for s in samples: s.runner.run_slurm(job_name = s.name+'_ray', **params)

########################## Link from Sisu to Taito ############################
old = "/homeappl/home/alice/"
new = "/wrk/bob/"
print "rsync -av --progress %s %s" % (proj.p.assembly_dir.path.replace(old, new), proj.p.assembly_dir)
for s in samples: print "rsync -av --progress %s %s" % (s.p.assembly_dir.path.replace(old, new), s.p.assembly_dir)

############################### Merged-Assembly ###############################
proj.merged.run(cpus=4)

################################## Mappings ###################################
params = dict(machines=1, cores=1, time='1-00:00:00', partition='serial',
              threads=6, mem_per_cpu=5300, constraint='hsw')
for s in samples: s.runner.run_slurm(steps=[{'mapper_51.run':{'cpus':6}}], job_name=s.name + "_eval_co_51_map", **params)
for s in samples: s.runner.run_slurm(steps=[{'mapper_61.run':{'cpus':6}}], job_name=s.name + "_eval_co_61_map", **params)
for s in samples: s.runner.run_slurm(steps=[{'mapper_71.run':{'cpus':6}}], job_name=s.name + "_eval_co_71_map", **params)
for s in samples: s.runner.run_slurm(steps=[{'mapper_81.run':{'cpus':6}}], job_name=s.name + "_eval_co_81_map", **params)
for s in samples: s.runner.run_slurm(steps=[{'mapper_merged.run':{'cpus':6}}], job_name=s.name + "_eval_m_map", **params)

params = dict(machines=1, cores=1, time='1-00:00:00', partition='serial',
              threads=6, mem_per_cpu=5300, constraint='hsw')
for s in samples: s.runner.run_slurm(steps=[{'mono_mapper.run':{'cpus':6}}], job_name=s.name + "_eval_mono_map",  **params)

################################# Binnings ####################################
params = dict(machines=1, cores=1, time='2-00:00:00', partition='longrun',
              threads=12, mem_per_cpu=5300, constraint='hsw')
proj.runner.run_slurm(steps=['assembly_51.results.binner.run'], job_name=proj.name+'_bin_51', **params)
proj.runner.run_slurm(steps=['assembly_61.results.binner.run'], job_name=proj.name+'_bin_61', **params)
proj.runner.run_slurm(steps=['assembly_71.results.binner.run'], job_name=proj.name+'_bin_71', **params)
proj.runner.run_slurm(steps=['assembly_81.results.binner.run'], job_name=proj.name+'_bin_81', **params)
proj.runner.run_slurm(steps=['merged.results.binner.run'],      job_name=proj.name+'_bin_04', **params)

################################## Plots ######################################
for s in tqdm(samples):
    print "Plots for sample '%s'" % s.name
    s.clean.fwd.graphs.length_dist.plot(x_log=True, y_log=True)
    s.clean.rev.graphs.length_dist.plot(x_log=True, y_log=True)
    s.assembly.results.contigs_fasta.graphs.length_dist.plot(x_log=True, y_log=True)
for a in proj.assemblies.values(): a.results.contigs_fasta.graphs.length_dist.plot(x_log=True, y_log=True)
proj.merged.results.contigs_fasta.graphs.length_dist.plot(x_log=True, y_log=True)

################################## Report #####################################
for s in tqdm(samples):
    print "Report on sample '%s'" % s.name
    s.report.generate()
