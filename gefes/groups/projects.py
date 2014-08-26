# Built-in modules #

# Internal modules #
from gefes.groups.aggregate import Aggregate
from gefes.groups.collection import Collection

# Third party modules #

###############################################################################
class Projects(Collection):
    """A collection of projects."""
    pass

###############################################################################
class Project(Aggregate):
    """A project containing several samples possibly spanning several runs."""

    def __repr__(self): return '<%s object "%s" with %i samples>' % \
                               (self.__class__.__name__, self.name, len(self))

    @property
    def long_name(self): return self.first.project_long_name

    def __init__(self, name, samples, projs_dir):
        # Attributes #
        self.name = name
        self.samples = samples
        # Dir #
        self.base_dir = projs_dir + self.name + '/'
        # Load #
        self.loaded = False