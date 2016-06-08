# -*- coding: utf-8 -*-

# Third party modules #
from ftputil import FTPHost
from tqdm import tqdm

###############################################################################
class SamplesToFTP(object):
    """Takes care of adding samples to a distant FTP server"""

    address = "ftp-private.ncbi.nlm.nih.gov"
    usrname = "subftp"
    passwrd = "w4pYB9VQ"
    dirctry = "/uploads/lucas.sinclair@me.com_Pm16RLfi"
    sub_dir = "PRJNA324769"

    def __init__(self, samples):
        # Save attributes #
        self.samples = samples

    def run(self, verbose=True):
        # Connect #
        if verbose: print "Connecting to FTP server..."
        self.ftp = FTPHost(self.address, self.usrname, self.passwrd)
        # Main loop #
        for sample in tqdm(self.samples):
            # Print #
            if verbose: print sample.short_name + ' (' + sample.name + ')'
            # Make names #
            self.base_name = sample.info['biosample'] + '_{}_reads.fastq.gz'
            # Test #
            assert sample.pair.fwd.count_bytes > 36
            assert sample.pair.rev.count_bytes > 36
            # Change to main directory #
            if verbose: print "Changing directories..."
            self.ftp.chdir(self.dirctry)
            # Make directory #
            if self.sub_dir not in self.ftp.listdir('.'):
                if verbose: print "Making directories..."
                self.ftp.mkdir(self.sub_dir)
            # Change to sub directory #
            self.ftp.chdir(self.sub_dir)
            # Upload #
            if verbose: print "Uploading forward..."
            self.ftp.upload(sample.pair.fwd, self.base_name.format("forward"))
            if verbose: print "Uploading reverse..."
            self.ftp.upload(sample.pair.rev, self.base_name.format("reverse"))
            # Return #
            self.ftp.close()
