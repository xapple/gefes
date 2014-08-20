# Built-in modules #
import os, xml.etree.ElementTree as etree

# Internal modules #
from gefes.groups.aggregate import Aggregate
from gefes.groups.collection import Collection
from plumbing.autopaths import AutoPaths

# Third party modules #

# Constants #
home = os.environ['HOME'] + '/'

###############################################################################
class Runs(Collection):
    """A collection of runs."""
    pass

###############################################################################
class Run(Aggregate):
    """An illumina run containing several pools."""

    def __repr__(self): return '<%s object number %i with %i pools>' % \
                               (self.__class__.__name__, self.num, len(self))

    def __init__(self, num, pools, out_dir):
        # Attributes #
        self.num = num
        self.name = "run%i" % num
        self.pools = pools
        # Dir #
        self.base_dir = out_dir + self.name + '/'
        self.p = AutoPaths(self.base_dir, self.all_paths)
        # Illumina report #
        self.xml_report_path = self.directory + "report.xml"
        # Auto exec #
        if os.path.isfile(self.xml_report_path): self.parse_report_xml()

    @property
    def label(self): return self.first.run_label

    @property
    def account(self): return self.first.account

    @property
    def directory(self):
        """The directory of the run"""
        return home + "proj/%s/INBOX/%s/" % (self.account, self.label)

    def parse_report_xml(self):
        tree = etree.parse(self.xml_report_path)
        root = tree.getroot()
        for sample in root.iter('Sample'):
            label = sample.items()[0][1]
            pool = [p for p in self if p.short_label == label][0]
            pool.report_stats['fwd'], pool.report_stats['rev'] = map(lambda x: dict(x.items()), sample.iter('Read'))