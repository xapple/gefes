
# Internal modules #
from gefes.common.autopaths import AutoPaths
from gefes.helper.contig import Contig
from gefes.helper.genecaller import GeneCaller

# Third party mods #
import sh

class BinAnnotater(object):
    """A bin is a object containing a multifasta of contigs that are ready to be annotated"""

    all_paths = """
    """

    def __init__(self,parent):
        self.parent = parent
        self.bini = parent
        
        
    

 
