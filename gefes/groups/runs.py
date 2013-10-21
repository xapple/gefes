# Built-in modules #
import xml.etree.ElementTree as etree

# Internal modules #
from aggregate import Collection, Aggregate
from gefes.common.autopaths import AutoPaths
import getpass

# Third party modules #

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
        # Parameters #
        self.label = self.first.run_label
        self.account = self.first.account
        # Dir #
        self.base_dir = out_dir + self.name + '/'
        self.p = AutoPaths(self.base_dir, self.all_paths)
        # Extra #
        if getpass.getuser() == 'inod':
            self.xml_report_path = "/proj/b2010008/INBOX/%s/report.xml" % (self.label)
        else:
            self.xml_report_path = "/proj/%s/INBOX/%s/report.xml" % (self.account, self.label)
        # Auto exec #
        self.parse_report_xml()

    def parse_report_xml(self):
        tree = etree.parse(self.xml_report_path)
        root = tree.getroot()
        for sample in root.iter('Sample'):
            label = sample.items()[0][1]
            pool = [p for p in self if p.short_label == label][0]
            pool.report_stats['fwd'], pool.report_stats['rev'] = map(lambda x: dict(x.items()), sample.iter('Read'))
