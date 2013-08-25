import hiseq

# Just one function for one pool #
pj = hiseq.projects['test']; p = pj[0]; p(steps=[{'assemble':{}}], threads=False)
# Just one pool via slurm #
hiseq.projects['test']['run000-pool01'].run_slurm()
hiseq.projects['humic']['run001-pool01'].run_slurm()

# Cut the first sample in half #
pj = hiseq.projects['humic']
p = pj[0]
p.fastq.cut_in_half(p.smaller)