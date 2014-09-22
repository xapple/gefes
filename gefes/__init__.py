b'This module needs Python 2.7.x'

# Special variables #
__version__ = '0.0.5'

# Built-in modules #
import os, sys, glob

# Get paths to module #
self = sys.modules[__name__]
module_dir = os.path.dirname(self.__file__)

# Dependencies #
from plumbing import dependencies
dependencies.check_setup(module_dir + 'setup.py')

# No need for an X display #
import matplotlib
matplotlib.use('Agg', warn=False)

# Internal modules #
from gefes.groups.samples import Sample
from gefes.groups.runs import Run, Runs
from gefes.groups.projects import Project, Projects
from plumbing.git import GitRepo

# Constants #
url = 'http://github.com/limno/gefes/'
home = os.environ['HOME'] + '/'

###############################################################################
# Output directory #
view_dir = home + 'GEFES/views/'

# Get paths to module #
repos_dir = os.path.abspath(module_dir + '/../') + '/'
git_repo = GitRepo(repos_dir)

# Load all samples #
json_paths = glob.glob(repos_dir + 'json/*/*.json')
samples = [Sample(j, view_dir + 'samples/') for j in json_paths]
samples.sort(key=lambda x: str(x))

# Compose into runs #
run_nums = sorted(list(set([s.run_num for s in samples])))
runs = [Run(num, [s for s in samples if s.run_num==num], view_dir + 'runs/') for num in run_nums]
runs = Runs(runs)
for s in samples: s.run = runs[s.run_num]

# Compose into projects #
proj_names = sorted(list(set([s.project_short_name for s in samples])))
projects = [Project(name, [s for s in samples if s.project_short_name==name], view_dir + 'projects/') for name in proj_names]
projects = Projects(projects)
for s in samples: s.project = projects[s.project_short_name]

# Make our favorite groups #
__import__("gefes.groups.favorites")