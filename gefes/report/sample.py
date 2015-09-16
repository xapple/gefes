# Futures #
from __future__ import division

# Built-in modules #
import os, json, socket
from collections import Counter, OrderedDict

# Internal modules #
import gefes

# First party modules #
from plumbing.autopaths import DirectoryPath, FilePath
from plumbing.common import split_thousands, pretty_now
from plumbing.cache import property_pickled
from pymarktex import Document, Template, HeaderTemplate, FooterTemplate
from pymarktex.figures import ScaledFigure, DualFigure

# Third party modules #
from tabulate import tabulate

# Constants #
ssh_header = "ssh://" + os.environ.get("FILESYSTEM_HOSTNAME", socket.getfqdn())
home = os.environ.get('HOME', '~') + '/'

###############################################################################
class SampleReport(Document):
    """A full report generated in PDF for every sample object."""

    def __repr__(self): return '<%s object on %s>' % (self.__class__.__name__, self.parent)

    def __init__(self, sample):
        # The parent #
        self.sample, self.parent = sample, sample
        # The output #
        self.base_dir    = self.sample.p.report_dir
        self.output_path = self.sample.p.report_pdf
        # Where should we cache stuff #
        self.cache_dir = DirectoryPath(self.base_dir + 'cached/')
        self.cache_dir.create(safe=True)

    def generate(self):
        # Dynamic templates #
        self.main = SampleTemplate(self)
        self.markdown = unicode(self.main)
        # Header and footer #
        self.header = HeaderTemplate()
        self.footer = FooterTemplate()
        # Render to latex #
        self.make_body()
        self.make_latex()
        self.make_pdf()

    uppmax_proj  = property(lambda self: self.sample.info.get('uppmax_project_id', 'b2014083'))
    export_base  = property(lambda self: 'GEFES/' + self.sample.project.name + '/' + self.sample.name + '.pdf')
    web_location = property(lambda self: FilePath(home + 'proj/' + self.uppmax_proj + '/webexport/' + self.export_base))
    url          = property(lambda self: "https://export.uppmax.uu.se/" + self.uppmax_proj + '/' + self.export_base)

###############################################################################
class SampleTemplate(Template):
    """All the parameters to be rendered in the markdown template"""
    delimiters = (u'{{', u'}}')

    def __init__(self, report):
        # Attributes #
        self.report, self.parent = report, report
        self.sample  = self.parent.sample
        self.project = self.sample.project
        self.cache_dir  = self.parent.cache_dir

    # General information #
    def sample_short_name(self):     return self.sample.name
    def sample_long_name(self):      return self.sample.long_name
    def project_short_name(self):    return self.project.name
    def project_long_name(self):     return self.project.long_name
    def project_other_samples(self): return len(self.project) - 1

    # Information dictionary #
    def json_url(self):
        json_location = os.path.relpath(self.project.json_path, start=gefes.git_repo)
        return gefes.url + ("tree/%s/" % gefes.git_repo.branch) + json_location
    def information(self):
        info = self.sample.info.copy()
        info.pop("contacts")
        info = json.dumps(info, indent=4, encoding='utf-8', ensure_ascii=False)
        info = info.strip("{}")
        return info

    # Process info #
    def project_url(self):       return gefes.url
    def project_version(self):   return gefes.__version__
    def git_hash(self):          return gefes.git_repo.hash
    def git_tag(self):           return gefes.git_repo.tag
    def git_branch(self):        return gefes.git_repo.branch
    def now(self):               return pretty_now()
    def results_directory(self): return ssh_header + self.sample.base_dir

    # Raw data #
    def fwd_size(self):        return             str(self.sample.pair.fwd.size)
    def rev_size(self):        return             str(self.sample.pair.rev.size)
    @property_pickled
    def fwd_count(self):       return split_thousands(self.sample.pair.fwd.count)
    @property_pickled
    def rev_count(self):       return split_thousands(self.sample.pair.rev.count)
    def illumina_report(self):
        if not self.sample.illumina_info.report.exists: return "<no report found>"
        else: return ssh_header + self.sample.illumina_info.report
    def raw_per_base_qual(self):
        params = [self.sample.pair.fwd.fastqc.results.per_base_qual,
                  self.sample.pair.rev.fastqc.results.per_base_qual]
        params += ["Forward",              "Reverse"]
        params += ["fwd_per_base_qual",    "rev_fwd_per_base_qual"]
        params += ["Raw per base quality", "raw_per_base_qual"]
        return str(DualFigure(*params))
    def raw_per_seq_qual(self):
        params = [self.sample.pair.fwd.fastqc.results.per_seq_qual,
                  self.sample.pair.rev.fastqc.results.per_seq_qual]
        params += ["Forward",                  "Reverse"]
        params += ["fwd_per_seq_qual",         "rev_per_seq_qual"]
        params += ["Raw per sequence quality", "raw_per_seq_qual"]
        return str(DualFigure(*params))

    # Preprocessing #
    def phred_threshold(self):   return self.sample.quality_checker.threshold
    def window_size(self):       return self.sample.quality_checker.window_size
    def length_threshold(self):  return self.sample.quality_checker.min_length
    @property_pickled
    def remaining_percent(self): return "%.2f%%" % (self.sample.quality_checker.results.ratio_kept * 100)
    @property_pickled
    def remaining_pairs(self):   return split_thousands(len(self.sample.quality_checker.dest))
    @property_pickled
    def remaining_singles(self): return split_thousands(len(self.sample.quality_checker.singletons))

    # Length distribution #
    def cleaned_len_dist(self):
        params = [self.sample.clean.fwd.length_dist.path,
                  self.sample.clean.rev.length_dist.path]
        params += ["Forward", "Reverse"]
        params += ["fwd_length_dist", "rev_length_dist"]
        params += ["Distribution of sequence lengths after quality control", "cleaned_len_dist"]
        return str(DualFigure(*params))
    def singletons_len_dist(self):
        caption = "Singletons length distribution"
        graph = self.sample.clean.singletons.length_dist
        label = "singletons_len_dist"
        return str(ScaledFigure(graph.path, caption, label))

    # FastQC graphs #
    def cleaned_per_base_qual(self):
        params = [self.sample.clean.fwd.fastqc.results.per_base_qual,
                  self.sample.clean.rev.fastqc.results.per_base_qual]
        params += ["Forward", "Reverse"]
        params += ["fwd_per_base_qual", "rev_per_base_qual"]
        params += ["Per base quality after quality control", "cleaned_per_base_qual"]
        return str(DualFigure(*params))
    def cleaned_per_seq_qual(self):
        params = [self.sample.clean.fwd.fastqc.results.per_seq_qual,
                  self.sample.clean.rev.fastqc.results.per_seq_qual]
        params += ["Forward", "Reverse"]
        params += ["fwd_per_seq_qual", "rev_per_seq_qual"]
        params += ["Per sequence quality after quality control", "cleaned_per_seq_qual"]
        return str(DualFigure(*params))

    # Rough taxonomic prediction #
    def kraken_version(self): return self.sample.kraken.long_name
    def kraken_domain_table(self):
        if not self.sample.kraken: return "<*Not computed yet*>"
        table = self.sample.kraken.results.at_domain_level
        table = OrderedDict((('Domain',table.keys()), ('Percentage',table.values())))
        table = tabulate(table, headers="keys", numalign="right", tablefmt="pipe")
        return table + "\n\n   : The Domain level breakdown predicted by Kraken."
    def kraken_phylum_table(self):
        if not self.sample.kraken: return "<*Not computed yet*>"
        table = self.sample.kraken.results.at_phylum_level
        table = OrderedDict((('Phylum',table.keys()), ('Percentage',table.values())))
        table = tabulate(table, headers="keys", numalign="right", tablefmt="pipe")
        return table + "\n\n   : The Phylum level distribution predicted by Kraken."
    def kraken_species_table(self):
        if not self.sample.kraken: return "<*Not computed yet*>"
        table = self.sample.kraken.results.at_species_level
        table = OrderedDict((('Species',table.keys()[0:10]), ('Percentage',table.values()[0:10])))
        table = tabulate(table, headers="keys", numalign="right", tablefmt="pipe")
        return table + "\n\n   : The 10 most common species predicted by Kraken."
    def kraken_summary_path(self): return ssh_header + self.sample.kraken.p.summary

    # Mono Assembly #
    def sample_assembler_version(self): return self.sample.assembly.long_name
    def sample_kmer_size(self):         return self.sample.assembly.kmer_size
    def sample_contig_cutoff(self):     return self.sample.assembly.length_cutoff
    def sample_count_contigs(self):
        if not self.sample.assembly: return 0
        return split_thousands(self.sample.assembly.results.contigs_fasta.count)
    def sample_contigs_len_dist(self):
        if not self.sample.assembly: return "<*Not computed yet*>"
        caption = "Mono-assembly length distribution"
        graph = self.sample.assembly.results.contigs_fasta.length_dist
        label = "sample_contigs_len_dist"
        return str(ScaledFigure(graph.path, caption, label))
    def sample_contigs_total_bp(self):
        if not self.sample.assembly: return 0
        return split_thousands(sum(self.sample.assembly.results.contigs_fasta.lengths))

    # Mono Mapping #
    def sample_mapper_version(self):   return self.sample.mono_mapper.long_name
    @property_pickled
    def sample_map_filter_count(self):
        if not self.sample.mono_mapper: return 0
        return split_thousands(self.sample.mono_mapper.results.filtered_count)
    @property_pickled
    def sample_did_map(self):
        if not self.sample.mono_mapper: return 0
        return "%.2f%%" % (self.sample.mono_mapper.results.fraction_mapped * 100)
    @property_pickled
    def sample_didnt_map(self):
        if not self.sample.mono_mapper: return 0
        return "%.2f%%" % (self.sample.mono_mapper.results.fraction_unmapped * 100)
    def sample_mean_coverage(self):
        if not self.sample.mono_mapper: return "<*Not computed yet*>"
        caption = "Mono-mapping mean coverage distribution"
        graph = self.sample.mono_mapper.results.mean_coverage_graph
        label = "sample_mean_coverage"
        return str(ScaledFigure(graph.path, caption, label))
    def samples_percent_covered(self):
        if not self.sample.mono_mapper: return "<*Not computed yet*>"
        caption = "Mono-mapping percent covered distribution"
        graph = self.sample.mono_mapper.results.percent_covered_graph
        label = "samples_percent_covered"
        return str(ScaledFigure(graph.path, caption, label))

    # Co Assembly #
    def count_contigs(self):    return split_thousands(self.sample.project.assembly.results.contigs_fasta.count)

    # Co-Mapping #
    def mapper_version(self):   return self.sample.mapper.long_name
    @property_pickled
    def map_filter_count(self):
        if not self.sample.mapper: return 0
        return split_thousands(self.sample.mapper.results.filtered_count)
    @property_pickled
    def did_map(self):
        if not self.sample.mapper: return 0
        return "%.2f%%" % (self.sample.mapper.results.fraction_mapped * 100)
    @property_pickled
    def didnt_map(self):
        if not self.sample.mapper: return 0
        return "%.2f%%" % (self.sample.mapper.results.fraction_unmapped * 100)
    def mean_coverage(self):
        if not self.sample.mapper: return "<*Not computed yet*>"
        caption = "Co-mapping mean coverage distribution"
        graph = self.sample.mapper.results.mean_coverage_graph
        label = "mean_coverage"
        return str(ScaledFigure(graph.path, caption, label))
    def percent_covered(self):
        if not self.sample.mapper: return "<*Not computed yet*>"
        caption = "Co-mapping percent covered distribution"
        graph = self.sample.mapper.results.percent_covered_graph
        label = "percent_covered"
        return str(ScaledFigure(graph.path, caption, label))

    # Protein calling (annotation) #
    def annotation_version(self):
        if not self.sample.contigs: return "<*Not computed yet*>"
        return self.sample.contigs[0].annotation.long_name
    @property_pickled
    def sample_count_proteins(self):
        if not self.sample.contigs: return "<*Not computed yet*>"
        if not all(c.annotation for c in self.sample.contigs): return "<*Not computed yet*>"
        total = sum(map(len,(c.annotation.results.functions for c in self.sample.contigs)))
        return split_thousands(total)
    @property_pickled
    def sample_functions_table(self):
        if not self.sample.contigs: return "<*Not computed yet*>"
        if not all(c.annotation for c in self.sample.contigs): return "<*Not computed yet*>"
        counts = Counter()
        for c in self.sample.contigs: counts.update(c.annotation.results.functions)
        table = OrderedDict(counts.most_common(20))
        table = {'Function': table.keys(), 'Counts': table.values()}
        table = tabulate(table, headers="keys", numalign="right", tablefmt="pipe")
        return table + "\n\n   : The 20 most common predicted functions in the predicted proteins of the mono-assembly."

