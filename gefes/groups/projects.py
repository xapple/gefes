# Built-in modules #
import re, os, glob

# Internal modules #
from gefes.groups.aggregates import Aggregate
from gefes.groups.collection import Collection
from gefes.groups.samples import Sample
from gefes.common import join_paired_filepaths

# First party modules #
from plumbing.autopaths import FilePath, AutoPaths
from plumbing.common import sort_string_by_pairs, natural_sort, load_json_path

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
    /info.json
    /samples/
    """

    def __repr__(self): return '<%s object "%s">' % \
                               (self.__class__.__name__, self.name)

    def __init__(self, json_path, project_dir):
        # Parse the json file describing the project #
        self.json_path = FilePath(json_path)
        self.info = load_json_path(json_path)
        # Required parameters #
        self.num       = self.info['project_num']
        self.name      = self.info['project_name']
        self.long_name = self.info['project_long_name']
        # Optional parameters #
        self.abstract  = self.info.get('abstract')
        self.run_name  = self.info.get('illumina_run_id')
        self.account   = self.info.get('uppmax_project_id')
        # Base directory #
        self.base_dir = project_dir + self.name + '/'
        # Delayed init #
        self.loaded = False

    def load(self):
        """A delayed kind of __init__ that is not called right away to avoid
        crowding the RAM of the python interpreter when you just import gefes"""
        # Load #
        self.loaded = True
        # Automatic paths #
        self.p = AutoPaths(self.base_dir, self.all_paths)
        # Make an alias to the json #
        self.json_path.link_to(self.p.info_json, safe=True)
        # Make all the samples object that this project possesses #
        if self.info.get('auto_parse_samples'):
            search_dir = self.info['samples_base_dir']
            search_dir = os.path.expanduser(search_dir)
            paths = glob.glob(search_dir + '*.fastq*')
            if not paths: paths = glob.glob(search_dir + '*.fasta*')
            if not paths: raise Exception("Found no FASTA or FASTQ path in %s" % search_dir)
            if all([re.search("_R[12]_", p) for p in paths]): pairs = join_paired_filepaths(paths)
            else:                                             pairs = sort_string_by_pairs(paths)
            self.samples = [Sample(self, p[0], p[1], num=i) for i,p in enumerate(pairs)]
            self.samples.sort(key=lambda x: natural_sort(x.name))
        else:
            self.samples = [Sample(self, info=info) for info in self.info['samples']]
        # The samples of a project are it's children in a way #
        self.children = self.samples
        # Call the mother function #
        return Aggregate.load(self)