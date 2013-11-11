# Built-in modules #
import os, gzip
from collections import Counter

# Internal modules #
from gefes.common.color import Color
from gefes.common.cache import property_cached
from gefes.common import isubsample
from gefes.common.autopaths import FilePath

# Third party modules #
import sh, shutil
from Bio import SeqIO

################################################################################
class FASTA(FilePath):
    """A single FASTA file somewhere in the filesystem"""
    extension = 'fasta'
    buffer_size = 1000

    def __len__(self): return self.count
    def __iter__(self): return self.parse()
    def __repr__(self): return '<%s object on "%s">' % (self.__class__.__name__, self.path)
    def __nonzero__(self): return os.path.getsize(self.path) != 0

    def __init__(self, path):
        self.path = path

    @property
    def first(self):
        self.open()
        seq = SeqIO.parse(self.handle, self.extension).next()
        self.close()
        return seq

    @property_cached
    def count(self):
        return int(sh.grep('-c', "^>", self.path, _ok_code=[0,1]))

    def open(self):
        if self.path.endswith('gz'): self.handle = gzip.open(self.path, 'r')
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

    def add_seq(self, seq):
        self.buffer.append(seq)
        self.buf_count += 1
        if self.buf_count % self.buffer_size == 0: self.flush()

    def flush(self):
        for seq in self.buffer: SeqIO.write(seq, self.handle, self.extension)
        self.buffer = []

    def write(self, reads):
        self.dir = os.path.dirname(self.path)
        if not os.path.exists(self.dir): os.makedirs(self.dir)
        self.handle = open(self.path, 'w')
        SeqIO.write(reads, self.handle, self.extension)
        self.handle.close()

    def copy(self, path):
        shutil.copy2(self.path, path)

    @property
    def lengths(self):
        return map(len, self.parse())

    @property_cached
    def lengths_counter(self):
        return Counter((len(s) for s in self.parse()))

    def subsample(self, down_to, new_path=None):
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

#-----------------------------------------------------------------------------#
class FASTQ(FASTA):
    """A single FASTQ file somewhere in the filesystem"""
    extension = 'fastq'

    @property_cached
    def count(self):
        return int(sh.grep('-c', "^+$", self.path, _ok_code=[0,1]))

    def to_fasta(self, path):
        with open(path, 'w') as handle:
            for r in self: SeqIO.write(r, handle, 'fasta')