# Futures #
from __future__ import division

# Built-in modules #
import os, json

# Internal modules #
from gefes.preprocess.quality import QualityChecker
from gefes.report.sample import SampleReport
from gefes.running.sample_runner import SampleRunner
from plumbing.autopaths import AutoPaths, FilePath
from fasta import PairedFASTQ
from fasta.fastqc import FastQC

# Third party modules #

# Constants #
home = os.environ['HOME'] + '/'

###############################################################################
class Sample(object):
    """Corresponds to an Illumina HiSeq MID. It's a bunch of paired sequences
    all coming from the same particular sample."""

    all_paths = """
    /info.json
    /clean/fwd.fastq
    /clean/rev.fastq
    /fastqc/fwd/
    /fastqc/rev/
    /mapping/
    /logs/
    /graphs/
    /report/report.pdf
    """

    def __repr__(self): return '<%s object "%s">' % (self.__class__.__name__, self.id_name)
    def __str__(self): return self.id_name
    def __iter__(self): return iter(self.children)
    def __len__(self): return self.count
    def __getitem__(self, key): return self.samples[key]

    def __init__(self, json_path, out_dir):
        # Output #
        self.out_dir = out_dir
        # Parse #
        self.json_path = FilePath(json_path)
        with open(json_path) as handle: self.info = json.load(handle)
        # Basic #
        self.account = self.info['uppmax_id']
        self.run_num = self.info['run_num']
        self.run_label = self.info['run_id']
        self.project_short_name = self.info['project']
        self.project_long_name = self.info['project_name']
        # Own attributes #
        self.num = self.info['pool_num']
        self.label = self.info['pool_id']
        self.short_label = self.label.split('_')[1]
        self.name = self.info['pool']
        self.long_name = self.info['pool_name']
        self.id_name = "run%03d-pool%02d" % (self.run_num, self.num)
        self.report_stats = {'fwd': {}, 'rev': {}}
        # Raw file pairs #
        fwd_path = home + "proj/%s/INBOX/%s/%s/%s" % (self.account, self.run_label, self.label, self.info['forward_reads'])
        rev_path = home + "proj/%s/INBOX/%s/%s/%s" % (self.account, self.run_label, self.label, self.info['reverse_reads'])
        self.pair = PairedFASTQ(fwd_path, rev_path)
        # Directory #
        self.base_dir = self.out_dir + self.id_name + '/'
        # Load #
        self.loaded = False

    def load(self):
        # Automatic paths #
        self.p = AutoPaths(self.base_dir, self.all_paths)
        # Runner #
        self.runner = SampleRunner(self)
        # Make an alias to the json #
        self.json_path.link_to(self.p.info_json, safe=True)
        # Change location of first FastQC #
        self.pair.fwd.fastqc = FastQC(self.pair.fwd, self.p.fastqc_fwd_dir)
        self.pair.rev.fastqc = FastQC(self.pair.rev, self.p.fastqc_rev_dir)
        # Cleaned pairs #
        self.clean = PairedFASTQ(self.p.fwd_clean, self.p.rev_clean)
        self.quality_checker = QualityChecker(self.pair, self.clean)
        self.singletons = self.quality_checker.singletons
        # Report #
        self.report = SampleReport(self)
        # Load #
        self.loaded = True
        # For convenience #
        return self

    @property
    def count(self):
        return self.pair.count

    def run_slurm(self, *args, **kwargs):
        return self.runner.run_slurm(*args, **kwargs)

    def clean_reads(self):
        self.cleaner.run()

    def map_reads(self):
        self.mapper.map()

    def make_plots(self):
        for graph in self.graphs: graph.plot()
