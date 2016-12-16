# Futures #
from __future__ import division

# Built-in modules #
import os, socket, shutil, inspect
from collections import OrderedDict

# Internal modules #
import gefes
from gefes.report import ReportTemplate

# First party modules #
from plumbing.autopaths import DirectoryPath, FilePath
from plumbing.common    import split_thousands, pretty_now
from plumbing.cache     import property_pickled
from pymarktex          import Document
from pymarktex.figures  import ScaledFigure

# Third party modules #
import pandas
from tabulate import tabulate

# Constants #
ssh_header = "ssh://" + os.environ.get("FILESYSTEM_HOSTNAME", socket.getfqdn())
home = os.environ.get('HOME', '~') + '/'

###############################################################################
class AssemblyReport(Document):
    """A full report generated in PDF for every assembly object."""

    def __repr__(self): return '<%s object on %s>' % (self.__class__.__name__, self.parent)

    def __init__(self, assembly):
        # The parent #
        self.assembly, self.parent = assembly, assembly
        self.aggregate = self.assembly.samples[0].project
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
        # Render to latex #
        self.make_body()
        self.make_latex()
        self.make_pdf()
        # Copy to reports directory #
        #shutil.copy(self.output_path, self.copy_base)
        # Return #
        return self.output_path

###############################################################################
class AssemblyTemplate(ReportTemplate):
    """All the parameters to be rendered in the markdown template"""
    delimiters = (u'{{', u'}}')

    def __init__(self, report):
        # Base parameters #
        self.report, self.parent = report, report
        self.assembly = self.parent.assembly
        self.cache_dir = self.parent.cache_dir
        # Other #
        self.aggregate = self.parent.aggregate

    # Assembly #
    def assembly_title(self):    return self.assembly.short_description
    def count_samples(self):     return len(self.assembly.samples)
    def assembler_version(self): return self.assembly.long_name
    def kmer_size(self):         return self.assembly.kmer_size
    def contig_cutoff(self):     return self.assembly.length_cutoff

    # General information #
    def aggregate_short_name(self): return self.aggregate.name
    def aggregate_long_name(self):  return self.aggregate.long_name

    # Process info #
    def results_directory(self): return ssh_header + self.aggregate.base_dir

    # Contigs #
    def count_contigs(self): return split_thousands(self.assembly.results.contigs_fasta.count)
    def contigs_len_hist(self):
        caption = "Assembly length distribution"
        graph = self.assembly.results.contigs_fasta.graphs.length_hist(x_scale='log', y_scale='log')
        label = "contigs_len_hist"
        return str(ScaledFigure(graph.path, caption, label))
    def contigs_total_bp(self): return split_thousands(self.assembly.results.total_bp)

    # Mapping #
    def mapping_version(self): return self.assembly.results.mappings.values()[0].long_name
    @property_pickled
    def mapping_table(self):
        info = OrderedDict((
             ('Name',      lambda s: "**" + s.name + "**"),
             ('Details',   lambda s: s.long_name),
             ('Reads',     lambda s: split_thousands(len(s.clean))),
             ('Did map',   lambda s: "%.1f%%" % (s.mappers[self.assembly].results.fraction_mapped * 100) if s.mappers[self.assembly] else "<*Not computed yet*>")
        ))
        row = lambda i: [f(self.assembly[i]) for f in info.values()]
        table = [[i+1] + row(i) for i in range(len(self.assembly))]
        table = tabulate(table, headers=info.keys(), numalign="right", tablefmt="pipe")
        return table + "\n\n   : Summary information for mapping of all samples."

    # Binning #
    def binning(self):
        if not self.assembly.results.binner: return False
        params = ('binning_version', 'count_bins', 'bins_contig_dist', 'bins_nucleotide_dist')
        return {p:getattr(self, p) for p in params}
    def binning_version(self): return self.assembly.results.binner.long_name
    def count_bins(self):      return split_thousands(len(self.assembly.results.binner.results))
    def bins_contig_dist(self):
        caption = "Bin number of contigs distribution"
        graph = self.assembly.results.binner.results.graphs.bins_contig_dist(x_scale='symlog')
        return str(ScaledFigure(graph.path, caption, inspect.stack()[0][3]))
    def bins_nucleotide_dist(self):
        caption = "Bin total nucleotide size distribution"
        graph = self.assembly.results.binner.results.graphs.bins_nucleotide_dist(x_scale='symlog')
        return str(ScaledFigure(graph.path, caption, inspect.stack()[0][3]))

    # Evaluation #
    def evaluation(self):
        if not self.assembly.results.binner: return False
        if not all(b.evaluation for b in self.assembly.results.binner.results.bins): return False
        params = ('bin_eval_version', 'bins_eval_graphs', 'bins_eval_markers_graph',
                  'bins_eval_marker_sets_graph', 'bins_eval_completeness_graph',
                  'bins_eval_contamination_graph', 'bins_eval_heterogeneity_graph',
                  'bins_eval_cch_graph', 'bins_quality_table', 'percent_mapped_to_good_bins',
                  'good_bins_count_contigs')
        return {p:getattr(self, p) for p in params}

    def bin_eval_version(self): return self.assembly.results.binner.results.bins[0].evaluation.long_name
    def bins_eval_graphs(self, name):
        graph = getattr(self.assembly.results.binner.results.eval_graphs, name)()
        return str(ScaledFigure(graph.path, "CheckMs '%s' metric" % name, "bins_eval_%s_graph" % name))

    def bins_eval_markers_graph(self):       return self.bins_eval_graphs('markers')
    def bins_eval_marker_sets_graph(self):   return self.bins_eval_graphs('marker_sets')
    def bins_eval_completeness_graph(self):  return self.bins_eval_graphs('completeness')
    def bins_eval_contamination_graph(self): return self.bins_eval_graphs('contamination')
    def bins_eval_heterogeneity_graph(self): return self.bins_eval_graphs('heterogeneity')
    def bins_eval_cch_graph(self):
        caption = "Contamination versus completeness with heterogeneity"
        graph = self.assembly.results.binner.results.eval_cch_graph()
        return str(ScaledFigure(graph.path, caption, "bins_eval_cch_graph"))

    # Quality bins #
    def bins_quality_table(self):
        """Sorted by completeness, report those that are >60%% complete and
        that have a contamination <10%. To this, add several other numbers:
        * The number of proteins found
        * The average coverage (across all samples)
        * The best taxonomic hit"""
        info = OrderedDict((('#',          lambda b: "**" + b.name + "**"),
                            ('Compl.',     lambda b: b.evaluation.results.statistics['completeness']),
                            ('Conta.',     lambda b: b.evaluation.results.statistics['contamination']),
                            ('Heter.',     lambda b: b.evaluation.results.statistics['heterogeneity']),
                            ('Prots.',     lambda b: len(b.faa)),
                            ('Avg. cov.',  lambda b: "%.2f" % b.average_coverage),
                            ('Assignment', lambda b: str(b.assignment.lowest_taxon))))
        good_bins = self.assembly.results.binner.results.good_bins
        frame = pandas.DataFrame(([f(b) for f in info.values()] for b in good_bins), columns=info.keys())
        frame = frame.set_index('#')
        frame = frame.sort_values(by="Compl.", ascending=False)
        table = tabulate(frame, headers='keys', numalign="right", tablefmt="pipe")
        return table + "\n\n   : Summary table for the best bins in this assembly."
    def percent_mapped_to_good_bins(self):
        frame        = self.assembly.results.mappings_per_sample
        good_bins    = self.assembly.results.binner.results.good_bins
        frame        = frame.loc[sum((b.contig_ids for b in good_bins), [])]
        total_reads  = sum(m.results.filtered_count for m in self.assembly.results.mappings.values())
        mapped_reads = frame.sum().sum()
        return "%.3f%%" % (100.0 * (float(mapped_reads) / float(total_reads)))
    def good_bins_count_contigs(self):
        return split_thousands(sum(map(len, self.assembly.results.binner.results.good_bins)))

    # Visualization #
    def visualization(self):
        if not self.assembly.results.binner.results.bins[-1].pfams: return False
        params = ('bin_genes_x_pfams', 'bin_bps_x_genes')
        return {p:getattr(self, p) for p in params}
    def bin_genes_x_pfams(self):
        caption = "Comparison of predicted number of genes against the number of predicted PFAMs"
        graph = self.assembly.results.binner.results.graphs.bin_genes_x_pfams()
        return str(ScaledFigure(graph.path, caption, inspect.stack()[0][3]))
    def bin_bps_x_genes(self):
        caption = "Comparison of cumulative length in base pairs against the number of predicted genes"
        graph = self.assembly.results.binner.results.graphs.bin_bps_x_genes()
        return str(ScaledFigure(graph.path, caption, inspect.stack()[0][3]))


