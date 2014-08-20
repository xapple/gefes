# Built-in modules #
from commands import getstatusoutput

# Internal modules #
from plumbing.autopaths import AutoPaths
from fasta import PairedFASTQ
from fasta import FASTQ

# Third party modules #

###############################################################################
class Pandaseq(object):
    """Takes care of running the pandaseq program which joins paired reads when they have overlap."""

    all_paths = """
    /lorem
    """

    def __repr__(self): return '<%s object of %s>' % (self.__class__.__name__, self.parent)
    def __len__(self): return len(self.pair)

    def __init__(self, cleaner):
        # Save parent #
        self.parent, self.cleaner = cleaner, cleaner
        self.pool = self.cleaner.pool
        # Auto paths #
        self.base_dir = self.parent.base_dir
        self.p = AutoPaths(self.base_dir, self.all_paths)
        # Files #
        pass

    def assemble(self):
        """Expects pandaseq 2.5"""
        command = 'pandaseq -f %s -r %s -u %s -F 1> %s 2> %s'
        command = command % (self.p.fwd_fastq, self.p.rev_fastq, self.unassembled.path, self.assembled.path, self.assembled.p.out)
        getstatusoutput(command)

    @property_cached
    def stats(self):
        result = {}
        result['raw'] = tail(self.p.out)
        if "pandaseq: error" in result['raw']: raise Exception("Pandaseq did not run properly")
        result['distrib'] = re.findall('STAT\tOVERLAPS\t(.+)$', result['raw'], re.M)
        result['distrib'] = map(int, result['distrib'][0].split())
        result['lengths'] = flatten([[i+1]*v for i,v in enumerate(result['distrib'])])
        result['noalign'] = int(re.findall('STAT\tNOALGN\t(.+)$', result['raw'], re.M)[0])
        result['lowqual'] = int(re.findall('STAT\tLOWQ\t(.+)$', result['raw'], re.M)[0])
        result['loss'] = 100 * sum(result['distrib'][100:]) / sum(result['distrib'])
        return result