# Built-in modules #
import sys
from collections import defaultdict

# Internal modules #
from gefes.binning import graphs
from gefes.binning.bin import Bin
from gefes.annotation.checkm import make_checkm_graphs, CheckmGraphCCH
from gefes.taxonomy.phylophlan import Phylophlan

# First party modules #
from plumbing.autopaths import AutoPaths
from plumbing.cache import property_cached
from plumbing.csv_tables import TSVTable

# Third party modules #
import pandas, sh
from tqdm import tqdm

###############################################################################
class Concoct(object):
    """Use CONCOCT at https://github.com/BinPro/CONCOCT
    to bin contigs together.
    Expects version 0.4.0.
    """

    short_name = 'concoct'
    long_name  = 'concoct v0.4.0'
    executable = 'concoct'
    dependencies = []

    all_paths = """
    /output/clustering_gt1000.csv
    /output/original_data_gt1000.csv
    /coverage.tsv
    /bins/
    /graphs/
    /taxonomy/
    """

    def __repr__(self): return '<%s object on %s>' % (self.__class__.__name__, self.assembly)
    def __nonzero__(self): return self.p.clustering.exists

    def __init__(self, samples, assembly, result_dir):
        # Save attributes #
        self.samples = samples
        self.assembly = assembly
        self.result_dir = result_dir
        # Auto paths #
        self.base_dir = self.result_dir + self.short_name + '/'
        self.p = AutoPaths(self.base_dir, self.all_paths)

    @property
    def coverage_matrix_tsv(self):
        """A dataframe where each row corresponds to a contig, and each column
        corresponds to a sample. The values are the average coverage for that contig
        in that sample."""
        tsv = TSVTable(self.p.coverage)
        if not tsv.exists:
            print "Computing coverage matrix for %i samples..." % len(self.samples); sys.stdout.flush()
            name_to_cov_mean = lambda s: self.assembly.results.mappings[s].results.coverage_mean
            frame = pandas.DataFrame({s.name: name_to_cov_mean(s.name) for s in self.samples})
            frame.to_csv(tsv.path, sep='\t', float_format='%.5g')
        return tsv

    @property_cached
    def coverage_matrix(self):
        return self.coverage_matrix_tsv.to_dataframe(index_col=0)

    def run(self):
        # Check samples loaded #
        for s in self.samples:
            if not s.loaded: s.load()
        # Check TSV is made #
        self.coverage_matrix_tsv
        # Run the pipeline #
        print "Launching CONCOCT..."; sys.stdout.flush()
        sh.concoct('--coverage_file',    self.coverage_matrix_tsv,
                   '--composition_file', self.assembly.results.contigs_fasta,
                   '-b',                 self.p.output_dir)
        # Remove the large and useless original data #
        self.p.original_data.remove()

    @property_cached
    def results(self):
        results = ConcoctResults(self)
        if not results: raise Exception("You can't access results from ConcoctResults before running the binning.")
        return results

###############################################################################
class ConcoctResults(object):

    def __nonzero__(self): return self.concoct.p.clustering.exists
    def __len__(self):     return len(self.bin_id_to_contig_ids)

    def __init__(self, concoct):
        self.concoct, self.parent = concoct, concoct
        self.p = concoct.p

    @property_cached
    def contig_id_to_bin_id(self):
        """Parse the raw result of CONOCT and return a dictionary with
        contig names as keys and bin names as values."""
        return dict(line.strip('\n').split(',') for line in self.concoct.p.clustering)

    @property_cached
    def bin_id_to_contig_ids(self):
        """The opposite of the above dictionary. Bin names as keys and
        lists of contig ids as values."""
        result = defaultdict(list)
        for c_id, b_id in self.contig_id_to_bin_id.items(): result[b_id].append(c_id)
        return result

    @property_cached
    def bins(self):
        """Return a list of all `Bin` objects."""
        return [Bin(self.concoct, c_ids, num=b_id) for b_id, c_ids in self.bin_id_to_contig_ids.items()]

    @property_cached
    def bin_id_to_bin(self):
        """A dictionary of bin ids to Bin objects."""
        return {b.name: b for b in self.bins}

    @property_cached
    def good_bins(self):
        """Return only the bins which are more than 60% complete and less than 10% contamination."""
        return [b for b in self.bins if b.evaluation.results.statistics['completeness']  > 60 and \
                                        b.evaluation.results.statistics['contamination'] < 10]

    @property_cached
    def good_contigs(self):
        """Return only the bins which are more than 60% complete and less than 10% contamination."""
        return [c for b in self.good_bins for c in b.contigs]

    @property_cached
    def taxonomy(self):
        """The results from running Phylophlan."""
        return Phylophlan(self.concoct, self.p.taxonomy_dir)

    @property_cached
    def graphs(self):
        class Graphs(object): pass
        result = Graphs()
        for graph in graphs.__all__:
            cls = getattr(graphs, graph)
            setattr(result, cls.short_name, cls(self))
        return result

    def run_all_bin_eval(self):
        """Run the evaluation procedure on all bins"""
        for b in tqdm(self.bins): b.evaluation.run()

    #-------------------------------------------------------------------------#
    @property_cached
    def eval_graphs(self):
        """An object with all the checkm graphs as attributes."""
        return make_checkm_graphs(self)

    @property_cached
    def eval_cch_graph(self):
        return CheckmGraphCCH(self)