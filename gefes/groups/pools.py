# Futures #
from __future__ import division

# Built-in modules #
import os, json

# Internal modules #
from gefes.common.autopaths import AutoPaths
from gefes.fasta.paired import PairedFASTQ
from gefes.fasta.single import FASTQ
from gefes.running.pool_runner import PoolRunner
from gefes.helper.cleaner import Cleaner
from gefes.graphs import pool_plots

# Third party modules #

###############################################################################
class Pool(object):
    """An illumina HiSeq MID is called here a 'pool'.
    It's a bunch of paired sequences."""

    all_paths = """
    /graphs/
    /clean/
    /info.json
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
        self.json_path = json_path
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
        self.short_name = self.info['pool']
        self.long_name = self.info['pool_name']
        self.id_name = "run%03d-pool%02d" % (self.run_num, self.num)
        self.report_stats = {'fwd': {}, 'rev': {}}
        # Automatic paths #
        self.base_dir = self.out_dir + self.id_name + '/'
        self.p = AutoPaths(self.base_dir, self.all_paths)
        # Make an alias to the json #
        if os.path.exists(self.p.info_json): os.remove(self.p.info_json)
        os.symlink(self.json_path, self.p.info_json)
        # Raw file pairs #
        self.fwd_path = "/proj/%s/INBOX/%s/%s/%s" % (self.account, self.run_label, self.label, self.info['forward_reads'])
        self.rev_path = "/proj/%s/INBOX/%s/%s/%s" % (self.account, self.run_label, self.label, self.info['reverse_reads'])
        # Convenience objects #
        self.fwd = FASTQ(self.fwd_path)
        self.rev = FASTQ(self.rev_path)
        self.pair = PairedFASTQ(self.fwd.path, self.rev.path)
        # Cleaning #
        self.cleaner = Cleaner(self)
        # All the plots #
        self.graphs = [getattr(pool_plots, cls_name)(self) for cls_name in pool_plots.__all__]
        # Runner #
        self.runner = PoolRunner(self)

    @property
    def count(self):
        return self.pair.count

    def run_slurm(self, *args, **kwargs):
        return self.runner.run_slurm(*args, **kwargs)

    def clean_reads(self):
        self.cleaner.clean()

    def make_plots(self):
        for graph in self.graphs: graph.plot()