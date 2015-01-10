# Built-in modules #
import os

# Internal modules #
from gefes.parsing.illumina      import IlluminaInfo
from gefes.preprocess.quality    import QualityChecker
from gefes.assemble.ray          import Ray
from gefes.map.bowtie            import Bowtie
from gefes.report.sample         import SampleReport
from gefes.running.sample_runner import SampleRunner
from gefes.annotation.prokka     import Prokka

# First party modules #
from plumbing.autopaths import AutoPaths, FilePath
from fasta              import PairedFASTA, PairedFASTQ
from fasta.fastqc       import FastQC

# Third party modules #

# Constants #
home = os.environ['HOME'] + '/'

###############################################################################
class Sample(object):
    """Consists of two FASTA or FASTQ files.
    It's a bunch of paired sequences all coming from the same particular lab sample.
    Might or not corresponds to an Illumina HiSeq MID. """

    all_paths = """
    /logs/
    /info.json
    /clean/fwd.fastq
    /clean/rev.fastq
    /clean/singletons.fastq
    /fastqc/fwd/
    /fastqc/rev/
    /assembly/
    /mapping/
    /graphs/
    /report/report.pdf
    """

    def __repr__(self): return '<%s object "%s">' % (self.__class__.__name__, self.name)
    def __iter__(self): return iter(self.children)
    def __len__(self): return self.count
    def __getitem__(self, key): return self.samples[key]

    def __init__(self, project, fwd_path=None, rev_path=None, out_dir=None, info=None, num=None, name=None):
        # Required parameters #
        self.project = project
        # Do we have the file paths already ? #
        if fwd_path is not None: self.fwd_path = FilePath(fwd_path)
        if rev_path is not None: self.rev_path = FilePath(rev_path)
        # Where should be put the output ? #
        if out_dir is None: self.out_dir = self.project.p.samples_dir
        else:               self.out_dir = out_dir
        # Let's inherit the information from the project #
        self.info = self.project.info.copy()
        if 'samples' in self.info: self.info.pop('samples')
        # Do we have extra information provided on this sample ? #
        if info is not None: self.info.update(info)
        # Extract fields from the extra information #
        self.long_name = self.info.get('samples_long_name')
        # Does the info contain the file location ? #
        self.raw_dir = self.info.get('samples_base_dir', '')
        if self.raw_dir and not self.raw_dir.endswith('/'): self.raw_dir += '/'
        self.raw_dir += self.info.get('sample_directory', '')
        if self.raw_dir and not self.raw_dir.endswith('/'): self.raw_dir += '/'
        if 'forward_reads' in self.info: self.fwd_path = FilePath(self.raw_dir + self.info.get('forward_reads'))
        if 'reverse_reads' in self.info: self.rev_path = FilePath(self.raw_dir + self.info.get('reverse_reads'))
        # Do we have a number for this sample ? #
        if num is None:  self.num = self.info.get('sample_num')
        else:            self.num = num
        # Do we have a name for this sample ? #
        if name is None: self.name = self.info.get('sample_name')
        else:            self.name = name
        # If we still don't have one, just use the file name #
        if self.name is None: self.name = self.fwd_path.short_prefix
        # Check that the files exist #
        if not self.fwd_path.exists: raise Exception("File '%s' does not exist" % self.fwd_path)
        if not self.rev_path.exists: raise Exception("File '%s' does not exist" % self.rev_path)
        # Is it a FASTA pair or a FASTQ pair ? #
        if "fastq" in self.fwd_path: self.pair = PairedFASTQ(self.fwd_path, self.rev_path)
        else:                        self.pair = PairedFASTA(self.fwd_path, self.rev_path)
        self.format = self.pair.format
        # Optional parameters #
        self.long_name = self.info.get('sample_long_name')
        # The directory where we will place all the data for this sample #
        self.base_dir = self.out_dir + self.name + '/'
        # Delayed init #
        self.loaded = False

    def load(self):
        """A delayed kind of __init__ that is not called right away to avoid
        crowding the RAM of the python interpreter when you just import gefes"""
        # Load #
        self.loaded = True
        # Check that the project is loaded #
        if not self.project.loaded: self.project.load()
        # Automatic paths #
        self.p = AutoPaths(self.base_dir, self.all_paths)
        # Maybe we have an Illumina report XML #
        self.illumina_info = IlluminaInfo(self)
        # Change location of first FastQC #
        if self.format == 'fastq':
            self.pair.fwd.fastqc = FastQC(self.pair.fwd, self.p.fastqc_fwd_dir)
            self.pair.rev.fastqc = FastQC(self.pair.rev, self.p.fastqc_rev_dir)
        # Cleaned pairs #
        if self.format == 'fastq':
            self.clean = PairedFASTQ(self.p.fwd_clean, self.p.rev_clean)
            self.quality_checker = QualityChecker(self.pair, self.clean)
            self.singletons = self.quality_checker.singletons
        # If it's a FASTA we can't clean it #
        if self.format == 'fasta': self.clean = self.pair
        # Map to the co-assembly #
        self.mapper = Bowtie(self, self.project.assembly, self.p.mapping_dir)
        # Assembly of this sample by itself #
        self.assembly = Ray([self], self.p.assembly_dir)
        # Annotate the contigs of the mono-assembly #
        self.annotation = Prokka(self, )
        # Runner #
        self.runner = SampleRunner(self)
        # Report #
        self.report = SampleReport(self)
        # For convenience #
        return self

    @property
    def count(self):
        return self.pair.count

    def run_slurm(self, *args, **kwargs):
        return self.runner.run_slurm(*args, **kwargs)

    def clean_reads(self):
        self.cleaner.run()

    def map_reads(self):
        self.mapper.map()

    def make_plots(self):
        for graph in self.graphs: graph.plot()
