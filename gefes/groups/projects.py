# Built-in modules #

# Internal modules #
from gefes.groups.aggregate import Aggregate
from gefes.groups.collection import Collection
from plumbing.autopaths import AutoPaths

# Third party modules #

###############################################################################
class Projects(Collection):
    """A collection of projects."""
    pass

###############################################################################
class Project(Aggregate):
    """A project containing several pools possibly spanning several runs."""

    def __repr__(self): return '<%s object "%s" with %i pools>' % \
                               (self.__class__.__name__, self.name, len(self))

    @property
    def long_name(self): return self.first.project_long_name

    def __init__(self, name, samples, projs_dir):
        # Attributes #
        self.name = name
        self.samples, self.children = samples, samples
        # Dir #
        self.base_dir = projs_dir + self.name + '/'

    def load(self):
        self.p = AutoPaths(self.base_dir, self.all_paths)
