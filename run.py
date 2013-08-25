# Just one function for one pool #
import hiseq; pj = hiseq.projects['test']; p = pj[0]; p(steps=[{'assemble':{}}], threads=False)
# Just one pool via slurm #
import hiseq; pj = hiseq.projects['test']; p = pj[0]; p.run_slurm()


# Cut the first sample in half #
ipython -i -c "import hiseq"
pj = hiseq.projects['humic']
p = pj[0]
p.fastq.cut_in_half(p.smaller)