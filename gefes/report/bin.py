# Futures #
from __future__ import division

# Built-in modules #
import os, socket, shutil, inspect

# Internal modules #
import gefes
from gefes.report import ReportTemplate

# First party modules #
from plumbing.autopaths import DirectoryPath, FilePath
from plumbing.cache     import property_cached
from pymarktex          import Document
from pymarktex.figures  import ScaledFigure

# Third party modules #

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
        self.assembly  = self.bin.assembly
        self.aggregate = self.assembly.samples[0].project
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
        # Copy to reports directory #
        shutil.copy(self.output_path, self.copy_base)
        # Return #
        return self.output_path

    @property_cached
    def copy_base(self):
        path = gefes.reports_dir + self.aggregate.name + '/bins/' + self.bin.name + '.pdf'
        return FilePath(path)

###############################################################################
class BinTemplate(ReportTemplate):
    """All the parameters to be rendered in the markdown template"""
    delimiters = (u'{{', u'}}')

    def __init__(self, report):
        # Attributes #
        self.report, self.parent = report, report
        self.bin        = self.parent.bin
        self.assembly   = self.parent.assembly
        self.project    = self.parent.aggregate
        self.cache_dir  = self.parent.cache_dir

    # General information #
    def bin_name(self):               return self.bin.name
    def bin_num(self):                return self.bin.num
    def binner_version(self):         return self.bin.binner.long_name
    def assembly_name(self):          return self.assembly.short_description
    def project_name(self):           return self.project.long_name
    def binner_other_bins(self):      return len(self.bin.binner.results.bins) - 1
    def bin_count_contigs(self):      return len(self.bin.contigs)
    def bin_count_prots(self):        return len(self.bin.faa)
    def good_bin_sentence(self):
        sentence = "It is %spart of the ``good bin'' group"
        return sentence % "" if self.bin.good else sentence % "not "

    # Process info #
    def results_directory(self): return ssh_header + self.bin.base_dir

    # Assignment #
    def assignment(self):
        return {'lowest_taxon': str(self.bin.assignment.lowest_taxon)} if self.bin.assignment else False

    # Visualization #
    def visualization(self):
        params = ('gc_x_totcov',)
        return {p:getattr(self, p) for p in params}
    def gc_x_totcov(self):
        caption = "Comparison of contig GC fraction against sum of average coverages"
        graph = self.bin.graphs.gc_x_totcov(rerun=True)
        return str(ScaledFigure(graph.path, caption, inspect.stack()[0][3]))

    # Annotation #
    def annotation(self): return False