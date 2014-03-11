
# Internal modules #
from gefes.common.autopaths import AutoPaths
from gefes.common.cache import property_cached
from gefes.fasta.single import FASTA
from gefes.helper.contig import Contig
from gefes.helper.genecaller import GeneCaller

class Bin(object):
    """A bin is a object containing a multifasta of contigs that are ready to be annotated"""

    all_paths = """
    /contigs.fasta
    /genes/
    """

    
    def __repr__(self): return '<%s object "%s" with %i Contigs>' % (self.__class__.__name__,self.name, len(self))
    def __iter__(self): return iter(self.contigs)
    def __len__(self): return len(self.contigs)
    def __getitem__(self, index): return self.contigs[index]
    
    def __init__(self,parent,c_list=[], name=None):
        self.parent = parent
        self.name = name
        self.contigs=c_list
        # Auto paths #
        self.base_dir = self.parent.p.bins + "/bin_" + name
        self.p = AutoPaths(self.base_dir, self.all_paths)
        self.caller=GeneCaller(self)
        
    def export(self):
        with FASTA(self.p.contigs) as ffile:
            for c in self.contigs: ffile.add_seq(c.record)

    def extend(self,contigs):
        self.contigs.extend(contigs)
                
    @classmethod
    def fromfolder(cls,parent, folder):
        #save parents and name
        name = folder.split("_")[1]
        bini=Bin(parent,[],name)
        # Auto paths #
        fasta = FASTA(bini.p.contigs)
        bini.extend([Contig(bini,f) for f in fasta])
        return bini
