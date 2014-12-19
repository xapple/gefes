# Built-in modules #

# Internal modules #
import xml.etree.ElementTree as etree

# First party modules #
from plumbing.autopaths import FilePath

# Third party modules #

###############################################################################
class IlluminaInfo(object):
    """Takes care of parsing Illumina XML reports."""

    def __init__(self, sample):
        # Base #
        self.sample = sample
        # Find the file #
        self.report = FilePath(self.sample.info.get("samples_base_dir") + "report.xml")

    @property
    def exists(self): return self.report.exists

    @property
    def stats(self):
        """Return stats_fwd and stats_rev """
        tree = etree.parse(self.report)
        root = tree.getroot()
        for sample in root.iter('Sample'):
            label = sample.items()[0][1]
            if label == self.sample.label: return map(lambda x: dict(x.items()), sample.iter('Read'))
        raise Exception("Could not find sample '%s' in the XML report '%s'" % (self.sample.name, self.report))

    @property
    def info(self):
        """The information for this sample"""
        if not self.exists: return {}
        return self.stats[0]