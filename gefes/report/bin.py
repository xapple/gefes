# Futures #
from __future__ import division

# Built-in modules #
import os, socket
from collections import Counter, OrderedDict

# Internal modules #
import gefes

# First party modules #
from plumbing.autopaths import DirectoryPath, FilePath
from plumbing.common import split_thousands, pretty_now
from plumbing.cache import property_pickled
from pymarktex import Document, Template
from pymarktex.figures import ScaledFigure, DualFigure

# Third party modules #
from tabulate import tabulate

# Constants #
ssh_header = "ssh://" + os.environ.get("FILESYSTEM_HOSTNAME", socket.getfqdn())
home = os.environ.get('HOME', '~') + '/'

###############################################################################
class BinReport(Document):
    """A full report generated in PDF for every bin object."""

    def __repr__(self): return '<%s object on %s>' % (self.__class__.__name__, self.parent)

    def __init__(self, bin):
        # The parent #
        self.bin, self.parent = bin, bin
        # The output #
        self.base_dir    = self.bin.p.report_dir
        self.output_path = self.bin.p.report_pdf
        # Where should we cache stuff #
        self.cache_dir = DirectoryPath(self.base_dir + 'cached/')
        self.cache_dir.create(safe=True)

    def generate(self):
        # Dynamic templates #
        self.main = BinTemplate(self)
        self.markdown = unicode(self.main)
        # Render to latex #
        self.make_body()
        self.make_latex()
        self.make_pdf()

###############################################################################
class BinTemplate(Template):
    """All the parameters to be rendered in the markdown template"""
    delimiters = (u'{{', u'}}')

    def __init__(self, report):
        # Attributes #
        self.report, self.parent = report, report
        self.bin        = self.parent.bin
        self.project    = self.sample.project
        self.cache_dir  = self.parent.cache_dir

    # General information #
    def bin_short_name(self):     return self.sample.name
    def bin_long_name(self):      return self.sample.long_name
    def assembly_short_name(self):    return self.project.name
    def assembly_long_name(self):     return self.project.long_name
    def assembly_other_samples(self): return len(self.project) - 1

    # Process info #
    def project_url(self):       return gefes.url
    def project_version(self):   return gefes.__version__
    def git_hash(self):          return gefes.git_repo.hash
    def git_tag(self):           return gefes.git_repo.tag
    def git_branch(self):        return gefes.git_repo.branch
    def now(self):               return pretty_now()
    def results_directory(self): return ssh_header + self.sample.base_dir