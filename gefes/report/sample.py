# Futures #
from __future__ import division

# Built-in modules #
import re, shutil

# Internal modules #
import gefes
from plumbing.common import split_thousands, pretty_now
from pymarktex import Document, Template, HeaderTemplate, FooterTemplate
from pymarktex.figures import DualFigure

# Third party modules #

###############################################################################
class SampleReport(Document):
    """A full report generated in PDF for every sample object."""

    def __repr__(self): return '<%s object on %s>' % (self.__class__.__name__, self.parent)

    def __init__(self, sample):
        self.sample, self.parent = sample, sample

    def generate(self):
        # Dynamic templates #
        self.md_template = SampleTemplate(self)
        self.header_template = HeaderTemplate()
        self.footer_template = FooterTemplate()
        # Render to latex #
        self.markdown = str(self.md_template)
        self.make_body()
        # Render to PDF #
        self.latex = str(self.header) + self.body + str(self.footer)
        self.make_pdf()

    def web_export(self):
        """Copy the report to the webexport directory where it can be viewed by anyone"""
        destination = "/proj/%s/webexport/GEFES/samples/run%03d_sample%02d.pdf"
        destination = destination % (self.sample.account, self.sample.run_num, self.sample.num)
        shutil.copy(self.p.pdf, destination)

    @property
    def url(self):
        link = "https://export.uppmax.uu.se/%s/ILLUMITAG/samples/run%03d_sample%02d.pdf"
        return link % (self.sample.account, self.sample.run_num, self.sample.num)

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
    def sample_long_name(self): return self.sample.short_name
    def project_short_name(self): return self.sample.project_short_name
    def project_long_name(self): return self.sample.project_long_name
    def project_other_samples(self): return len(self.project) - 1

    # JSON #
    def json_url(self):
        url = gefes.url + "/tree/master/json/"
        url += "/run%03d/run%03d-sample%03d.json"
        return url % (self.sample.run_num, self.sample.run_num, self.sample.num)
    def json_content(self):
        content = self.sample.json_path.read('utf-8')
        content = re.sub('\A(.+?)^    },$', '', content, flags=re.M|re.DOTALL)
        return content.strip('\n }')

    # Processing #
    def project_url(self): return gefes.url
    def project_version(self): return gefes.__version__
    def git_hash(self): return gefes.git_repo.hash
    def git_tag(self): return gefes.git_repo.tag
    def now(self): return pretty_now()
    def results_directory(self): return self.sample.base_dir

    # Raw data #
    def fwd_size(self): return str(self.sample.fwd.size)
    def fwd_count(self): return split_thousands(self.sample.fwd.count)
    def fwd_qual(self): return "%.2f" % self.sample.fwd.avg_quality
    def rev_size(self): return str(self.sample.rev.size)
    def rev_count(self): return split_thousands(self.sample.rev.count)
    def rev_qual(self): return "%.2f" % self.sample.rev.avg_quality
    def illumina_report(self): return self.sample.run.html_report_path
    def per_base_qual(self):
        params = [self.sample.fwd_fastqc.per_base_qual, self.sample.rev_fastqc.per_base_qual]
        params += ["Forward", "Reverse"]
        params += ["fwd_per_base_qual", "rev_fwd_per_base_qual"]
        params += ["Per base quality", "per_base_qual"]
        return str(DualFigure(*params))
    def per_seq_qual(self):
        params = [self.sample.fwd_fastqc.per_seq_qual, self.sample.rev_fastqc.per_seq_qual]
        params += ["Forward", "Reverse"]
        params += ["fwd_per_seq_qual", "rev_per_seq_qual"]
        params += ["Per sequence quality", "per_seq_qual"]
        return str(DualFigure(*params))

    # Filtering #
    def mismatches_allowed(self): return self.sample.assembled.primer_mismatches
    def primer_discard(self):
        before = self.sample.assembled
        after  = self.sample.assembled.good_primers.orig_reads
        return split_thousands(len(before) - len(after))
    def primer_left(self):
        return split_thousands(len(self.sample.assembled.good_primers.orig_reads))

    def n_base_discard(self):
        good = self.sample.assembled.good_primers
        return split_thousands(len(good.orig_reads) - len(good.n_filtered))
    def n_base_left(self):
        return split_thousands(len(self.sample.assembled.good_primers.n_filtered))

    def window_size(self): return self.sample.assembled.good_primers.qual_windowsize
    def window_threshold(self): return self.sample.assembled.good_primers.qual_threshold
    def window_discard(self):
        good = self.sample.assembled.good_primers
        return split_thousands(len(good.n_filtered) - len(good.qual_filtered))
    def window_left(self):
        return split_thousands(len(self.sample.assembled.good_primers.qual_filtered))

    def length_threshold(self): return self.sample.assembled.good_primers.min_length
    def length_discard(self):
        good = self.sample.assembled.good_primers
        return split_thousands(len(good.qual_filtered) - len(good.len_filtered))
    def length_left(self):
        return split_thousands(len(self.sample.assembled.good_primers.len_filtered))

    def percent_remaining(self):
        good = self.sample.assembled.good_primers
        return "%.1f%%" % ((len(good.len_filtered)/len(self.sample))*100)
