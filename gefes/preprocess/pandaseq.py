# Built-in modules #
import re

# Internal modules #
from plumbing.cache import property_cached
from plumbing.common import flatten, tail

# Third party modules #
from shell_command import shell_call

###############################################################################
class Pandaseq(object):
    """Takes care of running the pandaseq program which joins paired reads when they have overlap."""
    def __repr__(self): return "<%s object on '%s'>" % (self.__class__.__name__, self.source)

    def __init__(self, source, dest):
        self.source = source
        self.dest = dest

    def assemble(self):
        """We use shell_call because it exits with status 1.
        See https://github.com/neufeld/pandaseq/issues/40"""
        command = 'pandaseq27 -T 1 -f %s -r %s -u %s -F 1> %s 2> %s'
        command = command % (self.fwd, self.rev, self.unassembled.path, self.assembled.path, self.assembled.p.out)
        shell_call(command)

###############################################################################
class PandaseqResults(object):

    all_paths = """
    /cleaned_fwd.fastq
    /cleaned_rev.fastq
    /cleaned_single.fastq
    /report.txt
    """

    def __init__(self, base_dir):
        self.base_dir = base_dir

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