# Futures #
from __future__ import division

# Built-in modules #
import os, json, shutil, socket
from collections import Counter

# Internal modules #
import gefes
from plumbing.autopaths import FilePath
from plumbing.common import split_thousands, pretty_now
from pymarktex import Document, Template, HeaderTemplate, FooterTemplate
from pymarktex.figures import ScaledFigure, DualFigure

# Third party modules #
from tabulate import tabulate

# Constants #
ssh_header = "ssh://" + socket.getfqdn()

###############################################################################
class SampleReport(Document):
    """A full report generated in PDF for every sample object."""

    def __repr__(self): return '<%s object on %s>' % (self.__class__.__name__, self.parent)

    def __init__(self, sample):
        self.sample, self.parent = sample, sample
        self.output_path = self.sample.p.report_pdf

    def generate(self):
        # Dynamic templates #
        self.markdown = unicode(SampleTemplate(self))
        self.header = HeaderTemplate()
        self.footer = FooterTemplate()
        # Render to latex #
        self.make_body()
        self.make_latex()
        self.make_pdf()

    @property
    def location(self):
        loc = "GEFES/samples/run%03d_sample%02d.pdf"
        return loc % (self.sample.run_num, self.sample.num)

    def web_export(self):
        """Copy the report to the webexport directory where it can be viewed by anyone"""
        dest = FilePath(("/proj/%s/webexport/" + self.location) % self.sample.account)
        dest.make_directory()
        shutil.copy(self.output_path, dest)

    @property
    def url(self):
        return ("https://export.uppmax.uu.se/%s/" + self.location) % self.sample.account

###############################################################################
class SampleTemplate(Template):
    """All the parameters to be rendered in the markdown template"""
    delimiters = (u'{{', u'}}')

    def __init__(self, report):
        # Attributes #
        self.report, self.parent = report, report
        self.sample  = self.parent.sample
        self.project = self.sample.project

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
        info = json.dumps(info, sort_keys=True, indent=4, encoding='utf-8')
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
    def fwd_count(self):       return split_thousands(self.sample.pair.fwd.count)
    def rev_size(self):        return             str(self.sample.pair.rev.size)
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
    def remaining_percent(self): return "%.2f%%" % (self.sample.quality_checker.results.ratio_kept * 100)
    def remaining_pairs(self):   return split_thousands(len(self.sample.quality_checker.dest))
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
        graph = self.sample.clean.rev.length_dist
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

    # Mono Assembly #
    def sample_assembler_version(self): return self.sample.assembly.long_name
    def sample_kmer_size(self):         return self.sample.assembly.kmer_size
    def sample_contig_cutoff(self):     return self.sample.assembly.length_cutoff
    def sample_count_contigs(self):     return split_thousands(self.sample.assembly.results.contigs_fasta.count)
    def sample_contigs_len_dist(self):
        caption = "Mono-assembly length distribution"
        graph = self.sample.assembly.results.contigs_fasta.length_dist
        label = "sample_contigs_len_dist"
        return str(ScaledFigure(graph.path, caption, label))

    # Mono Mapping #
    def sample_mapper_version(self):   return self.sample.mono_mapper.long_name
    def sample_map_filter_count(self): return split_thousands(self.sample.mono_mapper.results.filtered_count)
    def sample_did_map(self):          return "%.2f%%" % (self.sample.mono_mapper.results.fraction_mapped * 100)
    def sample_didnt_map(self):        return "%.2f%%" % (self.sample.mono_mapper.results.fraction_unmapped * 100)
    def sample_mean_coverage(self):
        caption = "Mono-mapping mean coverage distribution"
        graph = self.sample.mono_mapper.results.mean_coverage_graph
        label = "sample_mean_coverage"
        return str(ScaledFigure(graph.path, caption, label))
    def samples_percent_covered(self):
        caption = "Mono-mapping percent covered distribution"
        graph = self.sample.mono_mapper.results.percent_covered_graph
        label = "samples_percent_covered"
        return str(ScaledFigure(graph.path, caption, label))

    # Protein calling (annotation) #
    def annotation_version(self): return self.sample.contigs[0].annotation.long_name
    def sample_functions(self):
        counts = Counter()
        for c in self.sample.contigs: counts.add(c.proteins)
        table = tabulate(counts.most_common(20), headers="keys", numalign="right", tablefmt="pipe")
        return table + "\n\n   : The 20 most abundant functions in the predicted proteins of the mono-assembly."