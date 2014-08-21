# Futures #
from __future__ import division

# Built-in modules #
import shutil


# Internal modules #
from gefes.common.autopaths import AutoPaths
from gefes.cleaning.sickle import Sickle
from gefes.fasta.single import FASTQ
from gefes.fasta.paired import PairedFASTQ

#from gefes.cleaning.quality import Quality 
#from gefes.cleaning.cutadapt import Cutadapt

# Third party modules #

###############################################################################
class Cleaner(object):
    """Takes care of cleaning the raw reads."""

    all_paths = """
    /fastqc/
    /cleaned_fwd.fastq
    /cleaned_rev.fastq
    """

    def __repr__(self): return '<%s object of %s>' % (self.__class__.__name__, self.parent)
    def __len__(self): return len(self.pair)

    def __init__(self, pool):
        # Save parent #
        self.parent, self.pool = pool, pool
        # Auto paths #
        self.base_dir = self.parent.p.clean_dir
        self.p = AutoPaths(self.base_dir, self.all_paths)
        self.sickle = Sickle(self)
#        self.quality = Quality(self)
#        self.cutadapt = Cutadapt(self)
#        shutil.copy(self.pool.fwd,self.p.fwd)
#        shutil.copy(self.pool.rev,self.p.rev)
        self.fwd=FASTQ(self.p.fwd)
        self.rev=FASTQ(self.p.rev)
        self.pair = PairedFASTQ(self.p.fwd, self.p.rev)

              
    @property
    def ratio_discarded(self):
        return 1 - (len(self.pair) / len(self.pool.pair))


    def run(self):
#        self.cutadapt.run()
        self.sickle.run()
#        self.quality.run()
