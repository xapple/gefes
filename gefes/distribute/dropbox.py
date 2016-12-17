# Built-in modules #

# Internal modules #
import gefes

# First party modules #

# Third party modules #
import sh

###############################################################################
class DropBoxSync(object):
    """Expects that your drop box is mounted at ~/Dropbox and then uses rsync.
    Don't forget to start the daemon:
        $ dropbox.py start
    """

    dbx_mount_point = gefes.home + "Dropbox/"

    def __init__(self, input_dir, output_dir):
        self.input_dir  = input_dir
        self.output_dir = output_dir

    def run(self):
        """Just rsync it"""
        print sh.rsync('-a', self.input_dir, gefes.home + 'Dropbox/' + self.output_dir)
