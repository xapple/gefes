# Built-in modules #
from collections import defaultdict

# Internal modules #
from gefes.binning.bin import Bin

# First party modules #
from plumbing.autopaths import AutoPaths
from plumbing.cache import property_cached
from plumbing.csv_tables import TSVTable

# Third party modules #
import pandas, sh

###############################################################################
class Concoct(object):
    """Use CONCOCT at https://github.com/BinPro/CONCOCT
    to bin contigs together.
    Expects version 0.4.0.
    """

    short_name = 'concoct'

    all_paths = """
    /output/clustering_gt1000.csv
    /output/original_data_gt1000.csv
    /coverage.tsv
    /bins/
    """

    def __repr__(self): return '<%s object on %s>' % (self.__class__.__name__, self.assembly)

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
            print "Computing coverage matrix for %i samples..." % len(self.samples)
            frame = pandas.DataFrame({s.name: s.mapper.results.coverage_mean for s in self.samples})
            frame.to_csv(tsv.path, sep='\t', float_format='%.5g')
        return tsv

    @property_cached
    def coverage_matrix(self):
        return self.coverage_matrix_tsv.to_dataframe(index_col=0)

    def run(self):
        # Run the pipeline #
        print "Launching CONCOCT..."
        sh.concoct('--coverage_file', self.coverage_matrix_tsv,
                   '--composition_file', self.assembly.results.contigs_fasta,
                   '-b', self.p.output_dir)
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
    def __init__(self, concoct):
        self.concoct = concoct

    @property_cached
    def contig_to_bin_id(self):
        """Parse the result of CONOCT and return a dictionary with
        contig names as keys and strings as values."""
        return dict(line.strip('\n').split(',') for line in self.concoct.p.clustering)

    @property_cached
    def bins(self):
        """Return a list of all `Bin` objects by making a dictionary with
        bin numbers as keys and a list of contig objects as values."""
        bins = defaultdict(list)
        for contig_name, bin_num in self.contig_to_bin_id.items():
            contig = [c for c in self.concoct.assembly.results.contigs if c.name == contig_name][0]
            bins[bin_num].append(contig)
        return [Bin(self.concoct, contigs, num=bin_num) for bin_num, contigs in bins.items()]