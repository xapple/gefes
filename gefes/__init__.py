b'This module needs Python 2.7.x'

# Special variables #
__version__ = '0.0.5'

# Built-in modules #
import os, sys, glob

# Get paths to module #
self = sys.modules[__name__]
module_dir = os.path.dirname(self.__file__)

# The module is a git repository #
from plumbing.git import GitRepo
repos_dir = os.path.abspath(module_dir + '/../') + '/'
git_repo = GitRepo(repos_dir)

# Dependencies #
from plumbing import dependencies
dependencies.check_setup(module_dir + 'setup.py')

# No need for an X display #
import matplotlib
matplotlib.use('Agg', warn=False)

# Internal modules #
from gefes.groups.projects import Project, Projects

# Constants #
url = 'http://github.com/limno/gefes/'
home = os.environ['HOME'] + '/'

###############################################################################
# Output directory #
view_dir = home + 'GEFES/views/'
project_dir = view_dir + "projects/"

# Load all projects by parsing project json files #
json_paths = glob.glob(repos_dir + 'json/*.json')
projects = [Project(j, project_dir) for j in json_paths]

# Convenience objects with indexing #
projects = Projects(projects)