# Built-in modules #
import os
from collections import OrderedDict

# Internal modules #
import gefes
from gefes.parsing.illumina      import IlluminaInfo
from gefes.preprocess.sliding    import SlidingWindow
from gefes.preprocess.sickle     import Sickle
from gefes.taxonomy.kraken       import Kraken
from gefes.assemble.ray          import Ray
from gefes.map.bowtie            import Bowtie
from gefes.report.sample         import SampleReport
from gefes.running.sample_runner import SampleRunner

# First party modules #
from plumbing.common    import load_json_path
from plumbing.autopaths import AutoPaths, FilePath, DirectoryPath
from plumbing.cache     import property_cached
from fasta              import PairedFASTA, PairedFASTQ
from fasta.fastqc       import FastQC

# Constants #
home = os.environ['HOME'] + '/'

###############################################################################
class Sample(object):
    """Consists of either two FASTA files or two FASTQ files.
    It's a bunch of paired sequences all coming from the same particular
    IRL lab sample. Might or might not correspond to an Illumina MID."""

    default_cleaner = "sickle"

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

    def __repr__(self):         return '<%s object "%s">' % (self.__class__.__name__, self.short_name)
    def __str__(self):          return self.short_name
    def __iter__(self):         return iter(self.children)
    def __len__(self):          return self.pair.count

    def __init__(self, json_path, raw_files_must_exist):
        # Attributes #
        self.json_path = FilePath(json_path)
        # Parse #
        self.info = load_json_path(self.json_path)
        self.info.pop('sentinel')
        # Own attributes #
        self.num                = self.info.get('sample_num')
        self.short_name         = self.info.get('sample_short_name')
        self.long_name          = self.info.get('sample_long_name')
        self.num                = int(self.info.get('sample_num'))
        # Project #
        self.project_short_name = self.info.get('project_short_name')
        self.project_long_name  = self.info.get('project_long_name')
        # Check the short name is only ASCII #
        assert all(ord(c) < 128 for c in self.short_name)
        assert self.short_name[0] not in "1234567890"
        # Automatic paths #
        self.base_dir  = gefes.samples_dir + self.info.get('organization') + '/'
        self.base_dir += self.project_short_name + '/' + self.short_name + '/'
        self.base_dir  = DirectoryPath(self.base_dir)
        self.p         = AutoPaths(self.base_dir, self.all_paths)
        # Make an alias to the json #
        self.p.info_json.link_from(self.json_path, safe=True)
        # Get the directory #
        prefix    = self.info.get('prefix',    '')
        directory = self.info.get('directory', '')
        suffix    = self.info.get('suffix',    '')
        # Get the file paths #
        self.fwd_path = FilePath(prefix + directory + suffix + self.info['fwd_filename'])
        self.rev_path = FilePath(prefix + directory + suffix + self.info['rev_filename'])
        # Is it a FASTA pair or a FASTQ pair ? #
        if "fastq" in self.fwd_path: self.pair = PairedFASTQ(self.fwd_path, self.rev_path)
        else:                        self.pair = PairedFASTA(self.fwd_path, self.rev_path)
        # Check that the files exist #
        if raw_files_must_exist:
            self.fwd_path.must_exist()
            self.rev_path.must_exist()
        # For speed let's update the sequence count cache if available #
        if self.info.get('forward_read_count') is not None:
            self.pair.fwd.count = self.info['forward_read_count']
        if self.info.get('reverse_read_count') is not None:
            self.pair.rev.count = self.info['reverse_read_count']
        # Change location of first FastQC, we don't want to modify the INBOX #
        if self.pair.format == 'fastq':
            self.pair.fwd.fastqc = FastQC(self.pair.fwd, self.p.fastqc_fwd_dir)
            self.pair.rev.fastqc = FastQC(self.pair.rev, self.p.fastqc_rev_dir)
        # Cleaned pairs: if it's a FASTA we can't clean it #
        if self.pair.format == 'fasta': self.clean = self.pair
        # Cleaned pairs: if it's a FASTQ #
        if self.pair.format == 'fastq': self.clean = PairedFASTQ(self.p.fwd_clean, self.p.rev_clean)

    #-------------------------------- Properties -----------------------------#
    @property_cached
    def illumina_info(self):
        """Maybe we have an Illumina report XML."""
        return IlluminaInfo(self)

    @property_cached
    def quality_checker(self):
        """With what are we going to preprocess the sequences and clean them ?"""
        assert self.pair.format == 'fastq'
        choices = {'window':   (SlidingWindow, (self.p.clean_dir, self.pair, self.clean)),
                   'sickle':   (Sickle,        (self.p.clean_dir, self.pair, self.clean))}
        cls, params = choices.get(self.default_cleaner)
        return cls(*params)

    @property_cached
    def mapper_51(self): return Bowtie(self, self.project.assembly_51, self.p.project_dir + "51/")
    @property_cached
    def mapper_61(self): return Bowtie(self, self.project.assembly_51, self.p.project_dir + "61/")
    @property_cached
    def mapper_71(self): return Bowtie(self, self.project.assembly_51, self.p.project_dir + "71/")
    @property_cached
    def mapper_81(self): return Bowtie(self, self.project.assembly_51, self.p.project_dir + "81/")
    @property_cached
    def mapper_merged(self): return Bowtie(self, self.project.merged, self.p.project_dir + "merged/")

    @property_cached
    def mappers(self):
        """A dictionary useful for trying different assemblies of different sizes.
        Keys are assembly objects, and values are corresponding mapper objects"""
        return OrderedDict(((self.project.assembly_51, self.mapper_51),
                            (self.project.assembly_61, self.mapper_61),
                            (self.project.assembly_71, self.mapper_71),
                            (self.project.assembly_81, self.mapper_81),
                            (self.project.merged, self.mapper_merged)))

    @property_cached
    def kraken(self):
        """Assembly of this sample by itself."""
        return Kraken(self.clean, self.p.kraken_dir)

    @property_cached
    def assembly(self):
        """Assembly of this sample by itself."""
        return Ray([self], self.p.assembly_dir)

    @property_cached
    def mono_mapper(self):
        """Map to the mono-assembly."""
        return Bowtie(self, self.assembly, self.p.mapping_mono_dir)

    @property_cached
    def runner(self):
        """The runner object."""
        return SampleRunner(self)

    @property_cached
    def report(self):
        """The PDF report."""
        return SampleReport(self)

    #-------------------------------- Shortcuts -----------------------------#
    @property
    def name(self):
        """Convenience shortcut. By default the short name."""
        return self.short_name

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

    @property
    def mapper(self):
        """Convenience shortcut. By default the merged assembly."""
        return self.mapper_merged
