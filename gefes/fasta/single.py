# Built-in modules #
import os, gzip
from collections import Counter, OrderedDict

# Internal modules #
from gefes.common import imean, isubsample
from gefes.common.color import Color
from gefes.common.cache import property_cached
from gefes.common.autopaths import FilePath, DirectoryPath
from gefes.common.tmpstuff import new_temp_dir, new_temp_path

# Third party modules #
import sh, shutil
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

################################################################################
class FASTA(FilePath):
    """A single FASTA file somewhere in the filesystem. You can read from it in
    several convenient ways. You can write to it in a automatically buffered way.
    There are several other things you can do with a FASTA file. Look at the class."""
    extension = 'fasta'
    buffer_size = 1000

    def __len__(self): return self.count
    def __iter__(self): return self.parse()
    def __repr__(self): return '<%s object on "%s">' % (self.__class__.__name__, self.path)
    def __nonzero__(self): return os.path.getsize(self.path) != 0
    def __contains__(self, other): return other in self.short_ids

    def __enter__(self): return self.create()
    def __exit__(self, exc_type, exc_value, traceback): self.close()

    def __getitem__(self, key):
        if isinstance(key, basestring): return self.sequences[key]
        elif isinstance(key, int): return self.sequences.items()[key]

    def __init__(self, path):
        self.path = path
        self.gziped = True if self.path.endswith('gz') else False

    @property
    def first(self):
        """Just the first sequence"""
        self.open()
        seq = SeqIO.parse(self.handle, self.extension).next()
        self.close()
        return seq

    @property_cached
    def count(self):
        """Should probably check file size instead of just caching once #TODO"""
        if self.gziped: return int(sh.zgrep('-c', "^+$", self.path, _ok_code=[0,1]))
        else: return int(sh.grep('-c', "^+$", self.path, _ok_code=[0,1]))

    @property_cached
    def short_ids(self):
        return frozenset([seq.description.split()[0] for seq in self])

    @property
    def lengths(self):
        return map(len, self.parse())

    @property_cached
    def lengths_counter(self):
        return Counter((len(s) for s in self.parse()))

    def open(self):
        if self.gziped: self.handle = gzip.open(self.path, 'r')
        else: self.handle = open(self.path, 'r')

    def close(self):
        if hasattr(self, 'buffer'): self.flush()
        self.handle.close()

    def parse(self):
        self.open()
        return SeqIO.parse(self.handle, self.extension)

    def create(self):
        self.buffer = []
        self.buf_count = 0
        self.dir = os.path.dirname(self.path)
        if not os.path.exists(self.dir): os.makedirs(self.dir)
        self.handle = open(self.path, 'w')
        return self

    def add_seq(self, seq):
        self.buffer.append(seq)
        self.buf_count += 1
        if self.buf_count % self.buffer_size == 0: self.flush()

    def add_str(self, seq, name=None, description=""):
        self.add_seq(SeqRecord(Seq(seq), id=name, description=description))

    def flush(self):
        for seq in self.buffer:
            SeqIO.write(seq, self.handle, self.extension)
        self.buffer = []

    def write(self, reads):
        self.dir = os.path.dirname(self.path)
        if not os.path.exists(self.dir): os.makedirs(self.dir)
        with open(self.path, 'w') as self.handle:
            SeqIO.write(reads, self.handle, self.extension)

    def get_id(self, id_num):
        """Extract one sequence from the file based on its ID. This is highly ineffective.
        Consider using the SQLite API instead."""
        for seq in self:
            if seq.id == id_num: return seq

    @property_cached
    def sequences(self):
        """Another way of easily retrieving sequences. Also highly ineffective.
        Consider using the SQLite API instead."""
        return OrderedDict(((seq.id,seq) for seq in self))

    def subsample(self, down_to=1, new_path=None):
        """Pick a number of sequences from the file randomly"""
        # Auto path #
        if not new_path: new_path = self.p.subsample
        # Check size #
        if down_to > len(self):
            message = "Can't subsample %s down to %i. Only down to %i."
            print Color.ylw + message % (self, down_to, len(self)) + Color.end
            self.copy(new_path)
            return
        # Make new file #
        self.subsampled = self.__class__(new_path)
        self.subsampled.create()
        # Do it #
        for seq in isubsample(self, down_to):
            self.subsampled.add_seqrecord(seq)
        # Clean up #
        self.subsampled.close()
        # Did it work #
        assert len(self.subsampled) == down_to

    def rename_with_num(self, prefix="", new_path=None, remove_desc=True):
        """Rename every sequence based on a prefix and a number"""
        # Temporary path #
        if new_path is None: numbered = self.__class__(new_temp_path())
        else: numbered = self.__class__(new_path)
        # Generator #
        def numbered_iterator():
            for i,read in enumerate(self):
                read.id = prefix + str(i)
                if remove_desc: read.description = ""
                yield read
        # Do it #
        numbered.write(numbered_iterator())
        numbered.close()
        # Replace it #
        if new_path is None:
            os.remove(self.path)
            shutil.move(numbered, self.path)

    def extract_length(self, lower_bound, upper_bound, new_path=None, cls=None):
        """Extract a certain length fraction and place them in a new file"""
        # Temporary path #
        cls = cls or self.__class__
        fraction = cls(new_temp_path()) if new_path is None else cls(new_path)
        # Generator #
        if lower_bound is None: lower_bound = 0
        def fraction_iterator():
            for read in self:
                if lower_bound <= len(read) <= upper_bound:
                    yield read
        # Do it #
        fraction.write(fraction_iterator())
        fraction.close()
        return fraction

#-----------------------------------------------------------------------------#
class FASTQ(FASTA):
    """A single FASTQ file somewhere in the filesystem"""
    extension = 'fastq'

    @property_cached
    def count(self):
        if self.gziped: return int(sh.zgrep('-c', "^+$", self.path, _ok_code=[0,1]))
        return int(sh.grep('-c', "^+$", self.path, _ok_code=[0,1]))

    def to_fasta(self, path):
        with open(path, 'w') as handle:
            for r in self: SeqIO.write(r, handle, 'fasta')
        return FASTA(path)

    def to_qual(self, path):
        with open(path, 'w') as handle:
            for r in self: SeqIO.write(r, handle, 'qual')
        return FilePath(path)

    @property_cached
    def avg_quality(self):
        mean = imean(s for r in self for s in r.letter_annotations["phred_quality"])
        self.close()
        return mean

    def fastqc(self, directory=None):
        # Default case #
        if directory is None:
            sh.fastqc(self.path, '-q')
            os.remove(self.prefix_path + '_fastqc.zip')
            return DirectoryPath(self.prefix.split('.')[0] + '_fastqc/')
        # Case directory #
        if directory is not None:
            if not isinstance(directory, DirectoryPath): directory = DirectoryPath(directory)
            if directory.exists: directory.remove()
            tmp_dir = new_temp_dir()
            sh.fastqc(self.path, '-q', '-o', tmp_dir)
            created_dir = tmp_dir + self.prefix.split('.')[0] + '_fastqc/'
            shutil.move(created_dir, directory)
            tmp_dir.remove()
            return directory