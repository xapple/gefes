# Built-in modules #
import os

# First party modules #
from fasta import FASTA
from plumbing.cache import property_cached

# Third party modules #
import ftputil

# Constants #
ftp_server = "ftp.ncbi.nlm.nih.gov"
home = os.environ['HOME'] + '/'
genome_dir = home + "/databases/genomes/"

###############################################################################
class RefGenome(object):
    """General reference genome object."""

    def __init__(self, source):
        self.source    = source
        self.file_name = os.path.basename(source)
        self.accession = os.path.splitext(self.file_name)[0]
        self.name      = os.path.basename(os.path.dirname(source))

    @property_cached
    def fasta(self):
        fasta = FASTA(genome_dir + self.name + "/genome.fasta")
        if fasta.exists: return fasta
        if not os.path.exists(genome_dir + self.name): os.makedirs(genome_dir + self.name)
        with ftputil.FTPHost("ftp.ncbi.nlm.nih.gov", "anonymous") as ftp: ftp.download(self.source, fasta)
        return fasta

###############################################################################
sources = ['/genomes/Bacteria/Escherichia_coli_K_12_substr__MG1655_uid57779/NC_000913.fna',
           '/genomes/Bacteria/Pseudomonas_putida_GB_1_uid58735/NC_010322.fna',
           '/genomes/Bacteria/Enterococcus_faecalis_62_uid159663/NC_017312.fna']

genomes = [RefGenome(s) for s in sources]