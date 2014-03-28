
#built-in modules#

# Internal modules #
from gefes.common.autopaths import AutoPaths
from gefes.common.slurm import nr_threads
from pandas import DataFrame

# Third party mods #
import sh
import os

class Phylotyper(object):
    """phylotyping bins with amphora of phylosift"""

    all_paths = """
    /phylosift/
    """

    sift = sh.Command(os.environ['PHYLOSIFT_PATH'])
    
    def __init__(self,parent):
        self.parent = parent
        self.bini = parent
        # Auto paths #
        self.base_dir = self.parent.p.phylotyping
        self.p = AutoPaths(self.base_dir, self.all_paths)

    def run(self):
        self.phylosift()

    def phylosift(self):
        self.sift("all", "-f", "--output=" + self.p.phylosift, self.parent.p.contigs)
                  

 
