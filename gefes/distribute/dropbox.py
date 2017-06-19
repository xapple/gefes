# Built-in modules #
import os, subprocess

# Internal modules #
import sifes

# First party modules #

# Third party modules #
import sh

# Constants #
home = os.environ.get('HOME', '~') + '/'

###############################################################################
class DropBoxRclone(object):
    """Expects that your dropbox is configured in rclone with the name "prod"
    Check with `$ rclone listremotes -l`."""

    def __init__(self, input_dir, output_dir):
        self.input_dir  = input_dir.rstrip('/')
        self.output_dir = output_dir.rstrip('/')

    @property
    def command(self):
        return 'sync', '--copy-links', "'" + self.input_dir + "'", "'" + 'prod:' + self.output_dir + "'"

    def run(self):
        """Just rclone it"""
        self.stdout = subprocess.check_call('rclone ' + ' '.join(self.command), shell=True)

###############################################################################
class DropBoxSync(object):
    """Expects that your dropbox is mounted at ~/Dropbox and then uses rsync locally.
    Don't forget to start the daemon:
        $ dropbox.py start
    """

    dbx_mount_point = home + "Dropbox/"

    def __init__(self, input_dir, output_dir):
        self.input_dir  = input_dir
        self.output_dir = output_dir

    def run(self):
        """Just rsync it and let the other script do the work."""
        assert "already running" in sh.Command("dropbox.py")("status")
        print sh.rsync('-a', self.input_dir, home + 'Dropbox/' + self.output_dir)
