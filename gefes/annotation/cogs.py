# Built-in modules #

# Internal modules #

# First party modules #
from plumbing.autopaths import AutoPaths

# Third party modules #

###############################################################################
class SingleCOGs(object):
    """Takes care of finding single copy COG genes in a set of contigs"""

    all_paths= """
    /lorem
    """

    def __init__(self, contigs, result_dir):
        # Save Attributes #
        self.contigs = contigs
        self.result_dir = result_dir
        # Auto paths #
        self.base_dir = self.result_dir + self.short_name + '/'
        self.p = AutoPaths(self.base_dir, self.all_paths)

    def run(self):
       pass