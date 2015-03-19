# Futures #
from __future__ import division

# Built-in modules #
import socket
from collections import OrderedDict

# Internal modules #
import gefes

# First party modules #
from plumbing.autopaths import DirectoryPath
from plumbing.common    import split_thousands, pretty_now
from plumbing.cache     import property_pickled
from pymarktex          import Document, Template, HeaderTemplate, FooterTemplate
from pymarktex.figures  import ScaledFigure

# Third party modules #
from tabulate import tabulate

# Constants #
ssh_header = "ssh://" + socket.getfqdn()

###############################################################################
class AssemblyReport(Document):
    """A full report generated in PDF for every assembly object."""

    def __repr__(self): return '<%s object on %s>' % (self.__class__.__name__, self.parent)

    def __init__(self, assembly):
        # The parent #
        self.assembly, self.parent = assembly, assembly
        # The output #
        self.base_dir    = self.assembly.p.report_dir
        self.output_path = self.assembly.p.report_pdf
        # Where should we cache stuff #
        self.cache_dir = DirectoryPath(self.base_dir + 'cached/')
        self.cache_dir.create(safe=True)

    def generate(self):
        # Dynamic templates #
        self.main = AssemblyTemplate(self)
        self.markdown = unicode(self.main)
        # Header and footer #
        self.header = HeaderTemplate()
        self.footer = FooterTemplate()
        # Render to latex #
        self.make_body()
        self.make_latex()
        self.make_pdf()

###############################################################################
class AssemblyTemplate(Template):
    """All the parameters to be rendered in the markdown template"""
    delimiters = (u'{{', u'}}')

    def __init__(self, report):
        # Base parameters #
        self.report, self.parent = report, report
        self.assembly = self.parent.assembly
        self.cache_dir = self.parent.cache_dir
        # Other #
        self.aggregate = self.assembly.samples[0].project

    # Assembly #
    def assembly_title(self):    return self.assembly.description
    def count_samples(self):     return len(self.assembly.samples)
    def assembler_version(self): return self.assembly.long_name
    def kmer_size(self):         return self.assembly.kmer_size
    def contig_cutoff(self):     return self.assembly.length_cutoff

    # General information #
    def aggregate_short_name(self): return self.aggregate.name
    def aggregate_long_name(self):  return self.aggregate.long_name

    # Process info #
    def project_url(self):       return gefes.url
    def project_version(self):   return gefes.__version__
    def git_hash(self):          return gefes.git_repo.hash
    def git_tag(self):           return gefes.git_repo.tag
    def git_branch(self):        return gefes.git_repo.branch
    def now(self):               return pretty_now()
    def results_directory(self): return ssh_header + self.aggregate.base_dir

    # Contigs #
    def count_contigs(self): return split_thousands(self.assembly.results.contigs_fasta.count)
    def contigs_len_dist(self):
        caption = "Assembly length distribution"
        graph = self.aggregate.assembly.results.contigs_fasta.graphs.length_dist.plot(x_log=True, y_log=True)
        label = "contigs_len_dist"
        return str(ScaledFigure(graph.path, caption, label))
    def contigs_total_bp(self): return split_thousands(sum(self.assembly.results.contigs_fasta.lengths))

    # Mapping #
    def mapping_version(self): return self.assembly.results.mappings.values()[0].long_name
    @property_pickled
    def mapping_table(self):
        info = OrderedDict((('Name',      lambda s: "**" + s.name + "**"),
                            ('Details',   lambda s: s.long_name),
                            ('Reads',     lambda s: split_thousands(len(s.clean))),
                            ('Did map',   lambda s: "%.3f%%" % (s.mappers[self.assembly].results.fraction_mapped * 100))))
        table = [[i+1] + [f(self.assembly[i]) for f in info.values()] for i in range(len(self.assembly))]
        table = tabulate(table, headers=info.keys(), numalign="right", tablefmt="pipe")
        return table + "\n\n   : Summary information for mapping of all samples."

    # Binning #
    def binning_version(self): return self.assembly.results.binner.long_name
    def count_bins(self):      return split_thousands(len(self.assembly.results.binner.results))
    def bins_contig_dist(self):
        caption = "Bin number of contigs distribution"
        graph = self.assembly.results.binner.results.graphs.bins_contig_dist(x_log=True)
        label = "bins_contig_dist"
        return str(ScaledFigure(graph.path, caption, label))
    def bins_nucleotide_dist(self):
        caption = "Bin total nucleotide size distribution"
        graph = self.assembly.results.binner.results.graphs.bins_nucleotide_dist(x_log=True)
        label = "bins_nucleotide_dist"
        return str(ScaledFigure(graph.path, caption, label))

    # Evaluation #
    def bin_eval_version(self): return self.assembly.results.binner.results.bins[0].evaluation.long_name
    def bins_genes_predicted_dist(self):
        caption = "Bin total nucleotide size distribution"
        graph = self.assembly.results.binner.results.graphs.bins_nucleotide_dist(x_log=True)
        label = "bins_nucleotide_dist"
        return str(ScaledFigure(graph.path, caption, label))
    def bins_eval_denisties_graph(self):
        caption = "Bin total nucleotide size distribution"
        graph = self.assembly.results.binner.results.graphs.bins_nucleotide_dist(x_log=True)
        label = "bins_nucleotide_dist"
        return str(ScaledFigure(graph.path, caption, label))
