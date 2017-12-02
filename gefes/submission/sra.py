# -*- coding: utf-8 -*-

"""
Example usage:

    import gefes
    from gefes.groups.lump import Lump
    projects = (bt, lb, kt)
    lump = Lump("ice_lump", projects)
    for b in tqdm(lump.good_bins): b.sra.upload_to_sra()
    lump.sra.write_sra_tsv()
"""

# Built-in modules #
import os, codecs

# Internal modules #
from plumbing.common    import Password, ascii
from plumbing.autopaths import AutoPaths

# Third party modules #
from ftputil import FTPHost

# Constants #
home = os.environ.get('HOME', '~') + '/'

# Old server #
ftp_server    = "ftp-private.ncbi.nih.gov"
ftp_login     = "sra"
ftp_password  = Password("SRA FTP password for user '%s':" % ftp_login)

# Lists #
header_bio = [
    '*sample_name',
    'description',
    'bioproject_id',
    'sample_title',
    '*organism',
    'host',
    'isolation_source',
    '*collection_date',
    '*env_biome',
    '*env_feature',
    '*env_material',
    '*geo_loc_name',
    '*lat_lon',
    'bin_number',
]

default_bio = {
    'bioproject_id'   : "PRJNA417934",
    'organism'        : "aquatic metagenome",
    'host'            : "not applicable",
    'isolation_source': "water",
    'depth'           : "surface",
    'env_biome'       : "lake",  # See https://trace.ddbj.nig.ac.jp/faq/biome-feature-material_e.html
    'env_feature'     : "lake",
    'env_material'    : "water",
}

header_sra =  [
    'bioproject_accession',
    'biosample_accession',
    'library_ID',
    'title',
    'library_strategy',
    'library_source',
    'library_selection',
    'library_layout',
    'platform',
    'instrument_model',
    'design_description',
    'genome_number',
    'filetype',
    'filename',
]

default_sra = {
    'bioproject_id'       : "PRJNA417044",
    'library_strategy'    : "AMPLICON",
    'library_source'      : "METAGENOMIC",
    'library_selection'   : "PCR",
    'library_layout'      : "Paired",
    'platform'            : "ILLUMINA",
    'instrument_model'    : "Illumina MiSeq",
    'forward_read_length' : "250",
    'reverse_read_length' : "250",
    'forward_filetype'    : "fastq",
    'reverse_filetype'    : "fastq",
}

###############################################################################
class BinSRA(object):
    """Every bin has an instance of this class, which enables us to
    access SRA required parameters for every MAG."""

    def __init__(self, bin):
        self.b = bin
        self.base_name = '%s_contigs.fastq'
        self.base_name = self.base_name % (self.b.long_name)
        self.first_sample = self.b.binner.samples[0]

    def upload_to_sra(self,
                      address = "ftp-private.ncbi.nlm.nih.gov",
                      usrname = "subftp",
                      passwrd = "w4pYB9VQ",
                      dirctry = "/uploads/lucas.sinclair@me.com_Pm16RLfi",
                      verbose = True):
        """They have an FTP site where you should drop the files first."""
        # Print #
        if verbose: print self.base_name
        # Connect #
        if verbose: print "Connecting..."
        self.ftp = FTPHost(address, usrname, passwrd)
        # Test #
        assert self.b.fasta > 36
        # Change directory #
        if verbose: print "Changing directories..."
        self.ftp.chdir(dirctry)
        # Make directory #
        if default_sra['bioproject_id'] not in self.ftp.listdir('.'):
            if verbose: print "Making directories..."
            self.ftp.mkdir(default_sra['bioproject_id'])
        self.ftp.chdir(default_sra['bioproject_id'])
        # Upload #
        if verbose: print "Uploading contigs..."
        self.ftp.upload(self.b.fasta, self.base_name)
        # Return #
        self.ftp.close()

    @property
    def biosample_line(self):
        """Will generate the corresponding BioSample entry (first TSV).
        Example usage:
            make_tsv.write_bio_tsv()
        You can add `from plumbing.common import gps_deg_to_float`
        """
        # sample_name #
        line = [self.b.long_name]
        # description #
        line += [self.b.project_long_name + ' bin number ' + self.b.num]
        # bioproject_id #
        line += [default_bio['bioproject_id']]
        # sample_title #
        line += ["Lake '%s' bin %s" % (self.b.long_name[-2:], self.b.num)]
        # organism #
        line += [default_bio["organism"]]
        # host #
        line += [default_bio["host"]]
        # isolation_source #
        line += [default_bio["isolation_source"]]
        # collection_date #
        line += [self.first_sample.info['date'].split(' ')[0]]
        # env_biome, env_feature, env_material #
        line += [default_bio["env_biome"]]
        line += [default_bio["env_feature"]]
        line += [default_bio["env_material"]]
        # geo_loc_name
        line += ["%s: %s" % (self.s.info['country'], ascii(self.s.info['location']))]
        # lat_lon
        coords = (float(self.first_sample.info['latitude'][0]),  # Latitude
                  float(self.first_sample.info['longitude'][0])) # Longitude
        line += ['{:7.6f} N {:7.6f} E'.format(*coords)]
        # bin_num
        line += [self.b.num]
        # Return #
        return line

    @property
    def sra_line(self):
        """Will generate the corresponding entry for SRA submission (second TSV).
        Example usage:
            make_tsv.write_sra_tsv()
        """
        # project accession
        bioproj = default_sra['bioproject_id']
        line    = [bioproj]
        # sample accession
        line += [self.first_sample.info['biosample']]
        # library id
        line += [self.b.long_name]
        # description
        desc    = "Lake '%s' metagenome assembled genome number %s"
        line   += [desc % (self.b.long_name[0:2], self.b.num)]
        # library_strategy
        line += [self.first_sample.info['library_strategy']]  # [default_sra['library_strategy']]
        line += [self.first_sample.info['library_source']]    # [default_sra['library_source']]
        line += [self.first_sample.info['library_selection']] # [default_sra['library_selection']]
        line += [self.first_sample.info['library_layout']]    #
        line += [self.first_sample.info['platform']]          # [default_sra['platform']]
        line += [self.first_sample.info['instrument_model']]  # [default_sra['instrument_model']]
        # design_description
        line += [self.first_sample.info['design_description']]
        # genome number
        line += [self.b.num]
        # general filetype
        line += ['fasta']
        # forward
        line += [self.base_name]
        # return
        return line

###############################################################################
class LumpSRA(object):
    """A class to generate the spreadsheets required by the NCBI/SRA
    for the submission of metagenome assembled genomes."""

    all_paths = """
    /bio_submission.tsv
    /sra_submission.tsv
    """

    def __init__(self, lump):
        """You give a lump as input."""
        # Base parameters #
        self.lump    = lump
        # Auto paths #
        self.base_dir = self.lump.p.sra_dir
        self.p        = AutoPaths(self.base_dir, self.all_paths)

    def write_bio_tsv(self, path=None):
        """Will write the TSV required by the NCBI for the creation of 'BioSample' objects
        (first TSV)."""
        # Content #
        header  = '\t'.join(header_bio) + '\n'
        content = '\n'.join('\t'.join(map(str, b.sra.sra_line)) for b in self.lump.good_bins)
        # Write it #
        if path is None: path = self.p.bio
        with codecs.open(path, 'w', encoding='utf-8') as handle:
            handle.write(header+content)
        # Message #
        print "Wrote TSV at '%s'" % path

    def write_sra_tsv(self, path=None):
        """Will write the appropriate TSV for the SRA submission in the cluster directory
        (second TSV). Sometimes you need to set the encoding to `windows-1252`."""
        # Content #
        header  = '\t'.join(header_sra) + '\n'
        content = '\n'.join('\t'.join(map(str, b.sra.sra_line)) for b in self.lump.good_bins)
        # Write it #
        if path is None: path = self.p.sra
        with codecs.open(path, 'w', encoding='utf-8') as handle:
            handle.write(header+content)
        # Message #
        print "Wrote TSV at '%s'" % path
