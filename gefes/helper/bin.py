# Internal modules #
from gefes.common.autopaths import AutoPaths
from gefes.fasta.single import FASTA
from gefes.helper.contig import Contig
from gefes.helper.genecaller import GeneCaller
from gefes.helper.bin_annotater import BinAnnotater
from gefes.running import Runner
from gefes.common.slurm import SLURMJob
from gefes.helper.phylotyper import Phylotyper
from gefes.helper.genome_builder import GenomeBuilder

################################################################################
class Bin(object):
    """A bin is a object containing a multifasta of contigs that are ready to be annotated"""

    all_paths = """
    /contigs.fasta
    /genes/
    /logs/
    /phylotyping/
    /rebuilt/
    """

    def __repr__(self): return '<%s object "%s" with %i Contigs>' % (self.__class__.__name__,self.name, len(self))
    def __iter__(self): return iter(self.contigs)
    def __len__(self): return len(self.contigs)
    def __getitem__(self, index): return self.contigs[index]

    def __init__(self,parent, c_list=[], name=None):
        self.parent = parent
        self.binner = parent
        self.name = name
        self.contigs=c_list
        # Auto paths #
        self.base_dir = self.parent.p.bins + "/bin_" + name
        self.p = AutoPaths(self.base_dir, self.all_paths)
        self.caller = GeneCaller(self)
        self.annotater = BinAnnotater(self)
        self.typer = Phylotyper(self)
        self.runner = BinRunner(self)
        self.genome_builder = GenomeBuilder(self)

    @classmethod
    def fromfolder(cls,parent, folder):
        #save parents and name
        name = folder.split("_")[1]
        bini=Bin(parent,[],name)
        # Auto paths #
        fasta = FASTA(bini.p.contigs)
        bini.extend([Contig(bini,f) for f in fasta])
        return bini

    def export(self):
        with FASTA(self.p.contigs) as ffile:
            for c in self.contigs: ffile.add_seq(c.record)

    def extend(self,contigs):
        self.contigs.extend(contigs)

    def calling(self):
        self.caller.run()

    def phylotyping(self):
        self.typer.run()

    def annotate(self):
        self.annotater.single_copy_cog_blast()

    def rebuild(self):
        self.genome_builder.pull_reads()
        self.genome_builder.assemble_genome()
        self.genome_builder.filter_assembly()
        

        
################################################################################

class BinRunner(Runner):
    """Will run stuff on a bin"""
    default_time = '12:00:00'

    default_steps = [
        {'phylotyping':{}},
        {'calling':    {}},
        {'annotate':   {}},
        {'rebuild':    {}}
    ]

    def __init__(self, parent):
        # Save parent #
        self.parent, self.bini = parent, parent
        self.binning = self.bini.parent
        self.project = self.binning.parent.parent

    def run_slurm(self, steps=None, **kwargs):
        # Make script #
        if not steps:
            steps = self.default_steps
        command = """steps = %s
                     binner = gefes.projects['%s'].binner
                     binner['%s'].load()
                     bini = [b for b in binner['%s'] if b.name=='%s'][0]
                     bini.runner(steps)""" % (steps,self.project.name, self.binning.name, self.binning.name,self.bini.name)
        # Test case #
        if 'test' in self.project.name:
            kwargs['time'] = '00:15:00'
            kwargs['qos'] = False
            kwargs['email'] = '/dev/null'
        # Send it #
        if 'time' not in kwargs: kwargs['time'] = self.default_time
        if 'email' not in kwargs: kwargs['email'] = None
        if 'job_name' not in kwargs: kwargs['job_name']  = "gefes_%s_%s_%s" % (self.project.name,self.binning.name, self.bini.name)
        self.slurm_job = SLURMJob(command, self.bini.p.logs_dir, **kwargs)
        self.slurm_job.launch()
