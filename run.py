# Just one function for one pool #
import hiseq; pj = hiseq.projects['test']; p = pj[0]; p(steps=[{'assemble':{}}], threads=False)
