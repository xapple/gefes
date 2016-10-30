# Built-in modules #

# Internal modules #
import gefes
from gefes.groups.aggregates import Aggregate
from gefes.groups.lump       import Lump

# Third party modules #

###############################################################################
class Projects(Lump):
    """A collection of projects."""
    def __init__(self, name="all_projects", *args, **kwargs):
        super(self.__class__,self).__init__(name, *args, **kwargs)

###############################################################################
class Project(Aggregate):
    """A project containing several samples. You can describe your
    projects in the JSON file placed in the repository root."""

    make_json_links = False

    all_paths = Aggregate.all_paths + """
    /info.json
    """

    def __init__(self, name, samples):
        """Please specify the name of the project and the samples it must contain."""
        # Base directory #
        out_dir = gefes.project_dir + self.organization + '/'
        # Super #
        super(self.__class__,self).__init__(name, samples, out_dir)

    @property
    def long_name(self): return self.first.project_long_name

    @property
    def organization(self): return self.firstinfo.get('organization', '')