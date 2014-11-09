# Built-in modules #
import json, glob

# Internal modules #
from gefes.groups.aggregates import Aggregate
from gefes.groups.collection import Collection
from gefes.groups.samples import Sample

# First party modules #
from plumbing.autopaths import FilePath, AutoPaths
from plumbing.common import sort_string_by_pairs

# Third party modules #

###############################################################################
class Projects(Collection):
    """A collection of projects."""
    pass

###############################################################################
class Project(Aggregate):
    """A project containing several samples. You can describe your
    projects in the a JSON placed in the repository."""

    all_paths = Aggregate.all_paths + """
    /samples/
    """

    def __repr__(self): return '<%s object "%s" with %i samples>' % \
                               (self.__class__.__name__, self.name, len(self))

    def __init__(self, json_path, project_dir):
        # Parse the json file describing the project #
        self.json_path = FilePath(json_path)
        with open(json_path) as handle: self.info = json.load(handle)
        # Required parameters #
        self.num       = self.info['project_num']
        self.name      = self.info['project_name']
        self.long_name = self.info['project_long_name']
        # Optional parameters #
        self.abstract  = self.info.get('abstract')
        # Base directory #
        self.base_dir = project_dir + self.name + '/'
        # Delayed init #
        self.loaded = False

    def load(self):
        """A delayed kind of __init__ that is not called right away to avoid
        crowding the RAM of the python interpreted when you just import gefes"""
        # Load #
        self.loaded = True
        # Automatic paths #
        self.p = AutoPaths(self.base_dir, self.all_paths)
        # Make all the samples object that this project posses #
        if self.info.get('auto_parse_samples'):
            search_dir = self.info['samples_base_dir']
            paths = glob.glob(search_dir + '*.fastq*')
            if not paths: paths = glob.glob(search_dir + '*.fasta*')
            if not paths: raise Exception("Found no FASTA or FASTQ path in %s" % search_dir)
            pairs = sort_string_by_pairs(paths)
            self.samples = [Sample(self, p[0], p[1], self.p.samples_dir, num=i) for i,p in enumerate(pairs)]
        else:
            message = "Ultimately you will be able to specify a list of file paths with extra metadata."
            raise NotImplementedError(message)
        # The samples of a project are it's children in a way #
        self.children = self.samples
        # Call the mother function #
        return Aggregate.load(self)