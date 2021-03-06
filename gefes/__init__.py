b'This module needs Python 2.7.x'

# Built-in modules #
import os, sys, socket

# Constants #
url         = 'http://xapple.github.io/gefes/'
repo_url    = 'http://github.com/xapple/gefes/'
__version__ = '1.0.2'
home        = os.environ.get('HOME', '~') + '/'
ssh_header  = "ssh://" + os.environ.get("FILESYSTEM_HOSTNAME", socket.getfqdn())

# No need for an X display #
import matplotlib
matplotlib.use('Agg', warn=False)

# Get paths to module #
self       = sys.modules[__name__]
module_dir = os.path.dirname(self.__file__)

# The module is maybe a git repository #
from plumbing.git import GitRepo
repos_dir = os.path.abspath(module_dir + '/../') + '/'
git_repo  = GitRepo(repos_dir, empty=True)

# Output directories #
view_dir    = home     + 'GEFES/views/'
project_dir = view_dir + 'projects/'
samples_dir = view_dir + 'samples/'
lumps_dir   = view_dir + 'lumps/'
reports_dir = home     + 'GEFES/reports/'
bundles_dir = home     + 'GEFES/bundles/'

# Internal modules #
from gefes.groups.projects import Projects

# The main objects, empty at first, call load() to populate them #
samples    = []
_projects  = []
projects   = Projects(aggregates=_projects)

# Expose functions #
from gefes.load import load