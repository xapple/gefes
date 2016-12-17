# Built-in modules #

# Internal modules #
import gefes
from gefes.groups.aggregates import Aggregate

# First party modules #
from plumbing.cache     import property_cached
from plumbing.autopaths import DirectoryPath

# Third party modules #

###############################################################################
class Bundle(Aggregate):
    """Regroup result files and reports from one or several projects for delivery."""

    all_paths = """
    /projects/
    /samples/
    /metadata/samples.xlsx
    """

    def __init__(self, name, samples, out_dir=None):
        # Directory #
        if out_dir is None: out_dir = gefes.bundles_dir
        # Super #
        super(self.__class__,self).__init__(name, samples, out_dir)
        # Figure out the projects within #
        proj_names = sorted(list(set([s.project_short_name for s in samples])))
        self.projects = [gefes.projects[p] for p in proj_names]

    def run(self):
        self.base_dir.remove()
        self.base_dir.create()
        for p in self.projects:
            # Directories #
            proj_dir = DirectoryPath(self.p.projects_dir + p.name)
            # Reports directory #
            reports_dir = DirectoryPath(proj_dir + 'reports')
            reports_dir.create(safe=True)
            # Reports for project #
            p.report.output_path.copy(reports_dir + 'project_report.pdf')
            # Reports for samples #
            smpl_rprts_dir = DirectoryPath(reports_dir + 'samples')
            smpl_rprts_dir.create(safe=True)
            for s in p: s.report.output_path.copy(smpl_rprts_dir + s.short_name + '.pdf')
            # Reports for assemblies #
            ably_rprts_dir = DirectoryPath(reports_dir + 'assemblies')
            ably_rprts_dir.create(safe=True)
            p.merged.results.report.output_path.copy(ably_rprts_dir + p.merged.short_name + '.pdf')
            # Data files #
            pass

    @property_cached
    def results(self):
        results = BundleResults(self)
        message = "You can't access results from a bundle before making the bundle."
        if not results: raise Exception(message)
        return results

###############################################################################
class BundleResults(object):

    def __nonzero__(self): return bool(self.p.x)

    def __init__(self, parent):
        self.parent   = parent
        self.base_dir = parent.base_dir
        self.p        = parent.p
