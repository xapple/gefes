# Futures #
from __future__ import division

# Built-in modules #
import socket
from collections import OrderedDict

# Internal modules #
import gefes

# First party modules #
from plumbing.autopaths import DirectoryPath
from plumbing.common import split_thousands, pretty_now
from plumbing.cache import property_pickled
from pymarktex import Document, Template, HeaderTemplate, FooterTemplate
from pymarktex.figures import ScaledFigure

# Third party modules #
from tabulate import tabulate

# Constants #
ssh_header = "ssh://" + socket.getfqdn()

###############################################################################
class AggregateReport(Document):
    """A full report generated in PDF for every aggregate object."""

    def __repr__(self): return '<%s object on %s>' % (self.__class__.__name__, self.parent)

    def __init__(self, aggregate):
        # The parent #
        self.aggregate, self.parent = aggregate, aggregate
        # The output #
        self.base_dir    = self.aggregate.p.report_dir
        self.output_path = self.aggregate.p.report_pdf
        # Where should we cache stuff #
        self.cache_dir = DirectoryPath(self.base_dir + 'cached/')
        self.cache_dir.create(safe=True)

    def generate(self):
        # Dynamic templates #
        self.main = AggregateTemplate(self)
        self.markdown = unicode(self.main)
        # Header and footer #
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
        self.cache_dir = self.parent.cache_dir

    # General information #
    def aggregate_short_name(self): return self.aggregate.name
    def aggregate_long_name(self):  return self.aggregate.long_name
    def count_samples(self):        return len(self.aggregate)

    # Samples #
    @property_pickled
    def sample_table(self):
        info = OrderedDict((
            ('Name',          lambda s: "**" + s.name + "**"),
            ('Details',       lambda s: s.long_name),
            ('Reads lost',    lambda s: "%.1f%%" % (100 - ((len(s.clean)/len(s.pair)) * 100))),
            ('Reads left',    lambda s: split_thousands(len(s.clean))),
            ('Mono mapped',   lambda s: "%.3f%%" % (s.mono_mapper.results.fraction_mapped * 100)),
            ('Co mapped',     lambda s: "%.3f%%" % (s.mapper.results.fraction_mapped * 100)),
        ))
        # The table #
        table = [[i+1] + [f(self.aggregate[i]) for f in info.values()] for i in range(len(self.aggregate))]
        # Make it as text #
        table = tabulate(table, headers=info.keys(), numalign="right", tablefmt="pipe")
        # Add caption #
        return table + "\n\n   : Summary information for all samples."

    # Process info #
    def project_url(self):       return gefes.url
    def project_version(self):   return gefes.__version__
    def git_hash(self):          return gefes.git_repo.hash
    def git_tag(self):           return gefes.git_repo.tag
    def git_branch(self):        return gefes.git_repo.branch
    def now(self):               return pretty_now()
    def results_directory(self): return ssh_header + self.aggregate.base_dir

    # Co Assembly #
    def assembler_version(self): return self.aggregate.assembly.long_name
    def kmer_size(self):         return self.aggregate.assembly.kmer_size
    def contig_cutoff(self):     return self.aggregate.assembly.length_cutoff
    def count_contigs(self):     return split_thousands(self.aggregate.assembly.results.contigs_fasta.count)
    def contigs_len_dist(self):
        caption = "Co-assembly length distribution"
        graph = self.aggregate.assembly.results.contigs_fasta.graphs.length_dist.plot(x_log=True, y_log=True)
        label = "contigs_len_dist"
        return str(ScaledFigure(graph.path, caption, label))

    # Binning #
    def binning_version(self): return self.aggregate.binner.long_name
    def count_bins(self):      return split_thousands(len(self.aggregate.binner.results))
    def bins_size_dist(self):
        caption = "Bin size distribution"
        graph = self.aggregate.binner.results.graphs.bins_size_dist.plot(x_log=True)
        label = "bins_size_dist"
        return str(ScaledFigure(graph.path, caption, label))
    def bins_nucleotide_dist(self):
        caption = "Bin size distribution"
        graph = self.aggregate.binner.results.graphs.bins_nucleotide_dist.plot(x_log=True)
        label = "bins_nucleotide_dist"
        return str(ScaledFigure(graph.path, caption, label))
