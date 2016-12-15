# Futures #
from __future__ import division

# Built-in modules #

# Internal modules #
from plumbing.common import flatter

# Third party modules #
import sh

###############################################################################
illumina_adapters = {
    'Adapter_1_Illumina'        : "ACACTCTTTCCCTACACGACGCTGTTCCATCT",
    'Adapter_4_Illumina'        : "CAAGCAGAAGACGGCATACGAGCTCTTCCGATCT",
    'Adapter_10_Illumina'       : "CAAGCAGAAGACGGCATACGAGATCGGTCTCGGCATTCCTGCTGAACCGCTCTTCCGATCT",
    'Adapter_12_Illumina'       : "CGGTCTCGGCATTCCTACTGAACCGCTCTTCCGATCT",
    'Adapter_15_Illumina'       : "ACAGGTTCAGAGTTCTACAGTCCGAC",
    'Adapter_18_Illumina'       : "AATGATACGGCGACCACCGACAGGTTCAGAGTTCTACAGTCCGA",
    'Adapter_19_Illumina'       : "CGACAGGTTCAGAGTTCTACAGTCCGACGATC",
    'Adapter_24_Illumina'       : "CCGACAGGTTCAGAGTTCTACAGTCCGACATG",
    'Adapter_26_Illumina'       : "TCGTATGCCGTCTTCTGCTTGT",
    'Adapter_32_Illumina'       : "TGGAATTCTCGGGTGCCAAGG",
    'Adapter_34_Illumina'       : "AGACGTGTGCTCTTCCGATC",
    'Adapter_36_Illumina'       : "AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT",
    'Adapter_37_Illumina'       : "GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT",
    'Adapter_53_Illumina'       : "GATCGTCGGACTGTAGAACTCTGAAC",
    'Adapter_75_TruSeq'         : "GATCGGAAGAGCACACGTCTGAACTCCAGTCACATCACGATCTCGTATGCCGTCTTCTGCTTG",
    'Adapter_76_TruSeq'         : "GATCGGAAGAGCACACGTCTGAACTCCAGTCACCGATGTATCTCGTATGCCGTCTTCTGCTTG",
    'Adapter_77_TruSeq'         : "GATCGGAAGAGCACACGTCTGAACTCCAGTCACTTAGGCATCTCGTATGCCGTCTTCTGCTTG",
    'Adapter_78_TruSeq'         : "GATCGGAAGAGCACACGTCTGAACTCCAGTCACTGACCAATCTCGTATGCCGTCTTCTGCTTG",
    'Adapter_79_TruSeq'         : "GATCGGAAGAGCACACGTCTGAACTCCAGTCACACAGTGATCTCGTATGCCGTCTTCTGCTTG",
    'Adapter_80_TruSeq'         : "GATCGGAAGAGCACACGTCTGAACTCCAGTCACGCCAATATCTCGTATGCCGTCTTCTGCTTG",
    'Adapter_81_TruSeq'         : "GATCGGAAGAGCACACGTCTGAACTCCAGTCACCAGATCATCTCGTATGCCGTCTTCTGCTTG",
    'Adapter_82_TruSeq'         : "GATCGGAAGAGCACACGTCTGAACTCCAGTCACACTTGAATCTCGTATGCCGTCTTCTGCTTG",
    'Adapter_83_TruSeq'         : "GATCGGAAGAGCACACGTCTGAACTCCAGTCACGATCAGATCTCGTATGCCGTCTTCTGCTTG",
    'Adapter_84_TruSeq'         : "GATCGGAAGAGCACACGTCTGAACTCCAGTCACTAGCTTATCTCGTATGCCGTCTTCTGCTTG",
    'Adapter_85_TruSeq'         : "GATCGGAAGAGCACACGTCTGAACTCCAGTCACGGCTACATCTCGTATGCCGTCTTCTGCTTG",
    'Adapter_86_TruSeq'         : "GATCGGAAGAGCACACGTCTGAACTCCAGTCACCTTGTAATCTCGTATGCCGTCTTCTGCTTG",
    'Adapter_87_TruSeq'         : "GATCGGAAGAGCACACGTCTGAACTCCAGTCACAGTCAACAATCTCGTATGCCGTCTTCTGCTTG",
    'Adapter_88_TruSeq'         : "GATCGGAAGAGCACACGTCTGAACTCCAGTCACAGTTCCGTATCTCGTATGCCGTCTTCTGCTTG",
    'Adapter_89_TruSeq'         : "GATCGGAAGAGCACACGTCTGAACTCCAGTCACATGTCAGAATCTCGTATGCCGTCTTCTGCTTG",
    'Adapter_90_TruSeq'         : "GATCGGAAGAGCACACGTCTGAACTCCAGTCACCCGTCCCGATCTCGTATGCCGTCTTCTGCTTG",
    'Adapter_91_TruSeq'         : "GATCGGAAGAGCACACGTCTGAACTCCAGTCACGTCCGCACATCTCGTATGCCGTCTTCTGCTTG",
    'Adapter_92_TruSeq'         : "GATCGGAAGAGCACACGTCTGAACTCCAGTCACGTGAAACGATCTCGTATGCCGTCTTCTGCTTG",
    'Adapter_93_TruSeq'         : "GATCGGAAGAGCACACGTCTGAACTCCAGTCACGTGGCCTTATCTCGTATGCCGTCTTCTGCTTG",
    'Adapter_94_TruSeq'         : "GATCGGAAGAGCACACGTCTGAACTCCAGTCACGTTTCGGAATCTCGTATGCCGTCTTCTGCTTG",
    'Adapter_95_TruSeq'         : "GATCGGAAGAGCACACGTCTGAACTCCAGTCACCGTACGTAATCTCGTATGCCGTCTTCTGCTTG",
    'Adapter_96_TruSeq'         : "GATCGGAAGAGCACACGTCTGAACTCCAGTCACGAGTGGATATCTCGTATGCCGTCTTCTGCTTG",
    'Adapter_97_TruSeq'         : "GATCGGAAGAGCACACGTCTGAACTCCAGTCACACTGATATATCTCGTATGCCGTCTTCTGCTTG",
    'Adapter_98_TruSeq'         : "GATCGGAAGAGCACACGTCTGAACTCCAGTCACATTCCTTTATCTCGTATGCCGTCTTCTGCTTG",
    'Adapter_99_Illumina'       : "GCCTTGGCACCCGAGAATTCCA",
    'Adapter_100_Illumina'      : "AATGATACGGCGACCACCGAGATCTACACGTTCAGAGTTCTACAGTCCGA",
    'Adapter_149_PCR_Primer1v2' : "AGATCGGAAGAGCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT",
    'Adapter_150_PCR_Primer1v3' : "GATCGGAAGAGCGTCGTGTAGGGGAAGAGGGTAGATCTCGGGGGGCGCCG",
}

###############################################################################
class Cutadapt(object):
    """Takes care of running the cutadapt program that will search for
    standard Illumina adapters that might be lingering in the sequences.
    https://github.com/marcelm/cutadapt
    There is a strange procedure for using paired-end files..."""

    def __repr__(self): return '<%s object of %s>' % (self.__class__.__name__, self.parent)
    def __len__(self): return len(self.pair)

    def __init__(self, source, dest):
        self.source = source
        self.dest = dest

    @property
    def adapter_params(self):
        return flatter([('-b', '%s=%s' % (k,v)) for k,v in illumina_adapters.items()])

    def run(self):
        sh.cutadapt(self.adapter_params + ['-O', 15, '-n', 2, '-o', self.p.cut_fwd, self.pool.fwd], _out=self.p.report_fwd)
        sh.cutadapt(self.adapter_params + ['-O', 15, '-n', 2, '-o', self.p.cut_rev, self.pool.rev], _out=self.p.report_rev)

###############################################################################
class CutadaptResults(object):

    all_paths = """
    /lorem
    """

    def __init__(self, base_dir):
        self.base_dir = base_dir