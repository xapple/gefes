# Built-in modules #

# Internal modules #
import sys, re
from collections import OrderedDict

# First party modules #
from plumbing.autopaths import AutoPaths
from plumbing.cache import property_cached
from plumbing.slurm import num_processors
from plumbing.graphs import Graph

# Third party modules #
import sh
from matplotlib import pyplot

###############################################################################
class Checkm(object):
    """Use CheckM at to evaluate a bin of contigs.
    Expects version v0.9.7.
    """

    short_name = 'checkm'
    long_name  = 'CheckM v0.9.7'
    executable = 'checkm'
    url        = 'https://github.com/Ecogenomics/CheckM'
    dependencies = ['hmmer', 'prodigal', 'pplacer']

    all_paths = """
    /contigs.fasta
    /stdout.txt
    /stderr.txt
    /output/
    """

    def __nonzero__(self): return bool(self.p.stdout)

    def __init__(self, bin, result_dir):
        # Save attributes #
        self.bin = bin
        # Auto paths #
        self.result_dir = result_dir
        self.p = AutoPaths(self.result_dir, self.all_paths)

    def run(self, cpus=None):
        # Variable threads #
        if cpus is None: cpus = num_processors
        # Link the bin's fasta file #
        self.bin.fasta.link_to(self.p.fasta)
        # Run the pipeline #
        print "Launching CheckM on bin '%s'..." % self.bin.name; sys.stdout.flush()
        sh.checkm('lineage_wf',
                  '-x', 'fasta',
                  '-t', cpus,
                  self.result_dir,
                  self.p.output_dir,
                  #'--tab_table', # See https://github.com/Ecogenomics/CheckM/issues/29
                  _out=self.p.stdout.path,
                  _err=self.p.stderr.path)
        # Check that it worked #
        assert 'unrecoverable error' not in self.p.stdout.contents

    @property_cached
    def results(self):
        results = CheckmResults(self)
        if not results: raise Exception("You can't access results from CheckM before running the algorithm.")
        return results

###############################################################################
class CheckmResults(object):

    def __nonzero__(self): return bool(self.checkm.p.stdout)

    def __init__(self, checkm):
        self.checkm = checkm

    @property_cached
    def statistics(self):
        """The various statistics produced by checkm in a dictionary. There is small technicality. The 'lineage' field values actually have a spaces in it, be careful when parsing. Use this rule: more than one space is necessary to split."""
        columns = OrderedDict((
            ("bin_id",      str),
            ("lineage",     str),
            ("genomes",     int),
            ("markers",     int),
            ("marker_sets", int),
            ("0", int), ("1", int), ("2", int), ("3", int), ("4", int), ("5+", int),
            ("completeness",  float),
            ("contamination", float),
            ("heterogeneity", float)))
        values = re.split(r'\s{2,}', list(self.checkm.p.stdout)[3])
        values = [v for v in values if v]
        return {k: columns[k](values[i]) for i,k in enumerate(columns)}

###############################################################################
def make_checkm_graphs(concot):
    """Will return an object with as attributes all the CheckM graphs.
    All graphs summarizing the results from the evaluation procedure.
    One graph for every statistic, later included in the assembly report."""
    # All the stats we want #
    names = ["genomes", "markers", "marker_sets",
             "completeness", "contamination", "heterogeneity"]
    # Make a dummy object #
    CheckmGraphs = type('CheckmGraphs', (), {})
    eval_graphs  = CheckmGraphs()
    # Main loop #
    for name in names:
        graph = CheckmSummaryGraph(concot, short_name=name)
        eval_graphs.__dict__[name] = graph
    return eval_graphs

class CheckmSummaryGraph(Graph):
    def plot(self, bins=250):
        counts = [b.evaluation.results.statistics.get(self.short_name) for b in self.parent.bins]
        fig = pyplot.figure()
        pyplot.hist(counts, bins=bins, color='gray')
        axes = pyplot.gca()
        axes.set_title("Distribution over all bins for the '%s' metric" % self.short_name)
        axes.set_xlabel(self.short_name)
        axes.set_ylabel('Number of bins with this much %s' % self.short_name)
        axes.xaxis.grid(False)
        self.save_plot(fig, axes, sep=('x'))
        pyplot.close(fig)
        return self