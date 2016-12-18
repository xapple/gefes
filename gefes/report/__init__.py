# Built-in modules #
import os, socket

# Internal modules #
import gefes

# First party modules #
from pymarktex import Template
from plumbing.common import pretty_now

# Constants #
ssh_header = "ssh://" + os.environ.get("FILESYSTEM_HOSTNAME", socket.getfqdn())

###############################################################################
class ReportTemplate(Template):
    """Things that are common to most reports in GEFES."""

    # Process info #
    def project_url(self):       return gefes.url
    def project_version(self):   return gefes.__version__
    def now(self):               return pretty_now()
    def git(self):
        if not gefes.git_repo: return False
        return {'git_hash'  : gefes.git_repo.hash,
                'git_tag'   : gefes.git_repo.tag,
                'git_branch': gefes.git_repo.branch}

    # Process info #
    def results_directory(self):
        return "ssh://cluster.sinclair.bio" + self.aggregate.base_dir
        return gefes.ssh_header + self.aggregate.base_dir
