# Futures #
from __future__ import division

# Built-in modules #
import re, shutil

# Internal modules #
import gefes
from plumbing.autopaths import FilePath
from plumbing.common import split_thousands, pretty_now
from pymarktex import Document, Template, HeaderTemplate, FooterTemplate
from pymarktex.figures import ScaledFigure, DualFigure

# Third party modules #

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
        self.sample = self.parent.sample
        self.project = self.sample.project

    # General information #
    def sample_short_name(self): return self.sample.short_name
    def sample_long_name(self): return self.sample.long_name
    def project_short_name(self): return self.sample.project_short_name
    def project_long_name(self): return self.sample.project_long_name
    def project_other_samples(self): return len(self.project) - 1

    # JSON #
    def json_url(self):
        url = gefes.url + "tree/%s/json/" % gefes.git_repo.branch
        url += "run%03d/run%03d-sample%03d.json"
        return url % (self.sample.run_num, self.sample.run_num, self.sample.num)
    def json_content(self):
        content = self.sample.json_path.read('utf-8')
        content = re.sub('\A(.+?)^    },$', '', content, flags=re.M|re.DOTALL)
        return content.strip('\n }')

    # Process info #
    def project_url(self): return gefes.url
    def project_version(self): return gefes.__version__
    def git_hash(self): return gefes.git_repo.hash
    def git_tag(self): return gefes.git_repo.tag
    def git_branch(self): return gefes.git_repo.branch
    def now(self): return pretty_now()
    def results_directory(self): return self.sample.base_dir

    # Raw data #
    def fwd_size(self):  return str(self.sample.pair.fwd.size)
    def fwd_count(self): return split_thousands(self.sample.pair.fwd.count)
    def rev_size(self):  return str(self.sample.pair.rev.size)
    def rev_count(self): return split_thousands(self.sample.pair.rev.count)
    def illumina_report(self): return self.sample.run.html_report_path
    def raw_per_base_qual(self):
        params = [self.sample.pair.fwd.fastqc.results.per_base_qual,
                  self.sample.pair.rev.fastqc.results.per_base_qual]
        params += ["Forward", "Reverse"]
        params += ["fwd_per_base_qual", "rev_fwd_per_base_qual"]
        params += ["Raw per base quality", "raw_per_base_qual"]
        return str(DualFigure(*params))
    def raw_per_seq_qual(self):
        params = [self.sample.pair.fwd.fastqc.results.per_seq_qual,
                  self.sample.pair.rev.fastqc.results.per_seq_qual]
        params += ["Forward", "Reverse"]
        params += ["fwd_per_seq_qual", "rev_per_seq_qual"]
        params += ["Raw per sequence quality", "raw_per_seq_qual"]
        return str(DualFigure(*params))

    # Preprocessing #
    def phred_threshold(self): return self.sample.quality_checker.threshold
    def window_size(self): return self.sample.quality_checker.window_size
    def length_threshold(self): return self.sample.quality_checker.min_length
    def remaining_percent(self): return "%.2f%%" % (self.sample.quality_checker.results.ratio_kept*100)
    def remaining_pairs(self): return split_thousands(len(self.sample.quality_checker.singletons))
    def remaining_singles(self): return split_thousands(len(self.sample.quality_checker.dest))

    # Length distribution #
    def cleaned_len_dist(self):
        params = [self.sample.clean.fwd.length_dist,
                  self.sample.clean.rev.length_dist]
        params += ["Forward", "Reverse"]
        params += ["fwd_length_dist", "rev_length_dist"]
        params += ["Distribution of sequence lengths after quality control", "cleaned_len_dist"]
        return str(DualFigure(*params))
    def singletons_len_dist(self):
        caption = "Singletons length distribution"
        path = self.sample.clean.rev.length_dist
        label = "singletons_len_dist"
        return str(ScaledFigure(path, caption, label))

    # Fastqc graphs #
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
