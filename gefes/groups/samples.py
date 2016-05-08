# Built-in modules #
import os, importlib
from collections import defaultdict, OrderedDict

# Internal modules #
from gefes.parsing.illumina      import IlluminaInfo
from gefes.preprocess.sliding    import SlidingWindow
from gefes.taxonomy.kraken       import Kraken
from gefes.assemble.ray          import Ray
from gefes.map.bowtie            import Bowtie
from gefes.report.sample         import SampleReport
from gefes.running.sample_runner import SampleRunner

# First party modules #
from plumbing.autopaths import AutoPaths, FilePath, DirectoryPath
from plumbing.cache     import property_cached
from fasta              import FASTQ, PairedFASTA, PairedFASTQ
from fasta.fastqc       import FastQC

# Third party modules #
from shell_command import shell_output
# Constants #
home = os.environ['HOME'] + '/'

###############################################################################
class Sample(object):
    """Consists of either two FASTA files or two FASTQ files.
    It's a bunch of paired sequences all coming from the same particular
    IRL lab sample. Might or might not correspond to an Illumina MID."""

    raw_files_must_exist = True

    all_paths = """
    /logs/
    /info.json
    /clean/fwd.fastq
    /clean/rev.fastq
    /clean/singletons.fastq
    /kraken/
    /fastqc/fwd/
    /fastqc/rev/
    /assembly/
    /mapping/project/
    /mapping/mono/
    /graphs/
    /report/report.pdf
    """

    def __repr__(self): return '<%s object "%s">' % (self.__class__.__name__, self.name)
    def __len__(self): return self.count

    def __init__(self, project, fwd_path=None, rev_path=None, out_dir=None, info=None, num=None, name=None):
        # Required parameters #
        self.project = project
        # Do we have the file paths already ? #
        if fwd_path is not None: self.fwd_path = FilePath(fwd_path)
        if rev_path is not None: self.rev_path = FilePath(rev_path)
        # Where should we put the output ? #
        if out_dir is None: self.out_dir = self.project.p.samples_dir
        else:               self.out_dir = out_dir
        # Let's inherit the information from the project #
        self.info = self.project.info.copy()
        if 'samples' in self.info: self.info.pop('samples')
        # Do we have extra information provided on this sample ? #
        if info is not None: self.info.update(info)
        # We should have the 'gefes_settings' key even if it is missing #
        if 'gefes_settings' not in self.info: self.info['gefes_settings'] = defaultdict(None)
        # Extract fields from the extra information #
        self.long_name = self.info.get('samples_long_name')
        # Does the info contain the file location ? #
        self.raw_dir = self.info.get('samples_base_dir', '')
        if self.raw_dir and not self.raw_dir.endswith('/'): self.raw_dir += '/'
        self.raw_dir += self.info.get('sample_directory', '')
        self.raw_dir = DirectoryPath(self.raw_dir)
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
        # Is it a FASTA pair or a FASTQ pair ? #
        if "fastq" in self.fwd_path: self.pair = PairedFASTQ(self.fwd_path, self.rev_path)
        else:                        self.pair = PairedFASTA(self.fwd_path, self.rev_path)
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
        # Check that the files exist #
        if self.raw_files_must_exist:
            self.fwd_path.must_exist()
            self.rev_path.must_exist()
        # For speed let's update the sequence count cache if available #
        if self.info.get('forward_read_count') is not None: self.pair.fwd.count = self.info['forward_read_count']
        if self.info.get('reverse_read_count') is not None: self.pair.rev.count = self.info['reverse_read_count']
        # Automatic paths #
        self.p = AutoPaths(self.base_dir, self.all_paths)
        # Maybe we have an Illumina report XML #
        self.illumina_info = IlluminaInfo(self)
        # Change location of first FastQC, we don't want to modify the INBOX #
        if self.pair.format == 'fastq':
            self.pair.fwd.fastqc = FastQC(self.pair.fwd, self.p.fastqc_fwd_dir)
            self.pair.rev.fastqc = FastQC(self.pair.rev, self.p.fastqc_rev_dir)
        # Cleaned pairs: if it's a FASTA we can't clean it #
        if self.pair.format == 'fasta': self.clean = self.pair
        # Cleaned pairs: if it's a FASTQ #
        if self.pair.format == 'fastq': self.clean = PairedFASTQ(self.p.fwd_clean, self.p.rev_clean)
        # Initial taxonomic predictions #
        self.kraken = Kraken(self.clean, self.p.kraken_dir)
        # Map to the co-assembly # TODO: rename the 71 mapper
        self.mapper = Bowtie(self, self.project.assembly, self.p.project_dir)
        # Map to different co-assemblies #
        self.mapper_51 = Bowtie(self, self.project.assembly_51, self.p.project_dir + "51/")
        self.mapper_61 = Bowtie(self, self.project.assembly_61, self.p.project_dir + "61/")
        self.mapper_71 = self.mapper
        self.mapper_81 = Bowtie(self, self.project.assembly_81, self.p.project_dir + "81/")
        self.mapper_merged = Bowtie(self, self.project.merged, self.p.project_dir + "merged/")
        # Assembly of this sample by itself #
        self.assembly = Ray([self], self.p.assembly_dir)
        # Map to the mono-assembly #
        self.mono_mapper = Bowtie(self, self.assembly, self.p.mapping_mono_dir)
        # Runner #
        self.runner = SampleRunner(self)
        # Report #
        self.report = SampleReport(self)
        # For convenience #
        return self

    #-------------------------------- Properties -----------------------------#
    @property_cached
    def quality_checker(self):
        """With what are we going to preprocess the sequences and clean them ?"""
        assert self.pair.format == 'fastq'
        info = self.info['gefes_settings'].get('quality_checker')
        if info is None: return SlidingWindow(self.p.clean_dir, self.pair, self.clean)
        module = importlib.import_module(info['object']['source'])
        obj    = getattr(module, info['object']['name'])
        return obj(self.p.clean_dir, self.pair, self.clean)

    @property
    def mappers(self):
        """A dictionary useful for trying different assemblies of different sizes.
        Keys are assembly objects, and values are corresponding mapper objects"""
        return OrderedDict(((self.project.assembly_51, self.mapper_51),
                            (self.project.assembly_61, self.mapper_61),
                            (self.project.assembly_71, self.mapper_71),
                            (self.project.assembly_81, self.mapper_81),
                            (self.project.merged, self.mapper_merged)))

    #------------------------------ Special cases ---------------------------#
    def merge_lanes(self, remove_orig=False):
        """We got a few runs that had several lanes in the same sample directory.
        We want to `cat` all these to files called fwd.fastq.gz and rev.fastq.gz"""
        # Find all lanes #
        fwd_match = lambda f: f.endswith('R1_001.fastq.gz')
        rev_match = lambda f: f.endswith('R2_001.fastq.gz')
        fwd_files = [FASTQ(f) for f in self.raw_dir.flat_files if fwd_match(f)]
        rev_files = [FASTQ(f) for f in self.raw_dir.flat_files if rev_match(f)]
        for x,y in zip(fwd_files, rev_files): print "Combining these files:", x.prefix, 'with', y.prefix
        shell_output("zcat %s |gzip > %s" % (' '.join(fwd_files), self.pair.fwd))
        shell_output("zcat %s |gzip > %s" % (' '.join(rev_files), self.pair.rev))
        # Check #
        assert self.pair.fwd.first.id == self.pair.rev.first.id
        # Remove the original #
        if remove_orig:
            for f in fwd_files + rev_files: f.remove()

    #-------------------------------- Shortcuts -----------------------------#
    @property
    def count(self):
        """Convenience shortcut. The number of sequences of the raw pair."""
        return self.pair.count

    @property
    def singletons(self):
        """Convenience shortcut. The singletons of the quality check."""
        return self.quality_checker.singletons

    @property
    def contigs(self):
        """Convenience shortcut. The contigs of the mono-assembly."""
        if not self.assembly: return []
        return self.assembly.results.contigs