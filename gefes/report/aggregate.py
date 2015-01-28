# Futures #
from __future__ import division

# Built-in modules #
import re, shutil

# Internal modules #
import gefes
from plumbing.common import split_thousands, pretty_now
from pymarktex import Document, Template, HeaderTemplate, FooterTemplate
from pymarktex.figures import ScaledFigure, DualFigure

# Third party modules #

###############################################################################
class AggregateReport(Document):
    """A full report generated in PDF for every aggregate object."""

    def __repr__(self): return '<%s object on %s>' % (self.__class__.__name__, self.parent)

    def __init__(self, sample):
        # The parent #
        self.aggregate, self.parent = sample, sample
        # The output #
        self.base_dir    = self.sample.p.report_dir
        self.output_path = self.sample.p.report_pdf

    def generate(self):
        # Dynamic templates #
        self.markdown = unicode(AggregateTemplate(self))
        self.header = HeaderTemplate()
        self.footer = FooterTemplate()
        # Render to latex #
        self.make_body()
        self.make_latex()
        self.make_pdf()

###############################################################################
class AggregateTemplate(Template):
    """All the parameters to be rendered in the markdown template"""
    delimiters = (u'{{', u'}}')

    def __init__(self, report):
        # Attributes #
        self.report, self.parent = report, report
        self.aggregate = self.parent.aggregate

    # General information #
    def sample_short_name(self):     return self.sample.name
    def sample_long_name(self):      return self.sample.long_name
    def project_short_name(self):    return self.project.name
    def project_long_name(self):     return self.project.long_name
    def project_other_samples(self): return len(self.project) - 1
