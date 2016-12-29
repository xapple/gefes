# Futures #
from __future__ import division

# Built-in modules #
from collections import defaultdict

# Internal modules #
import inspect, os

# First party modules #
from plumbing.cache     import property_cached
from plumbing.autopaths import AutoPaths

# Third party modules #
import pandas

# Current directory #
filename    = inspect.getframeinfo(inspect.currentframe()).filename
current_dir = os.path.dirname(os.path.abspath(filename)) + '/'

###############################################################################
class HitProfile(object):
    """The idea is to build a profile from protein matches (such as Pfam or Tigrfam) that shows the
    different metabolic processes occurring along a profile of samples.
    The samples can be either a time series or geospatial series.

    #---------------------#
    1. Input: contigs_filtered.fasta, BAM files, and pfam/hits.hmmout file.

    `samples_x_contigs`
    2. Matrix of average coverage per contigs.
    ROW: contigs - COL: samples - VAL: Average coverage.

    `contigs_x_pfams`
    3. Matrix of HMM search results.
    ROW: pfams   - COL: contigs - VAL: Count observations.
    -> For making matrix 3, only take the top hit of each predicted protein.

    `samples_x_pfams`
    4: Multiply M3*M2: Samples x Pfams.
    ROW: samples  - COL: pfams    - VAL: Sums of average coverages.
    -> This will 'merge' columns which have have the same Pfam by addition.

    `normalization_vector`
    5. Select the columns that correspond to the single copy Pfams. Calculate median by row (sample).
    -> Normalize by row using this number.

    `norm_samples_x_pfams`
    6. Divide M4 with S5.
    -> Final outputed matrix for depth profile.

    Pitfalls: Coverage variation over contigs (check with IVG browser).
    Pitfalls: PFAM subfamilies (the extra number is the version of the PFAM model)
              can be thrown away and results merged by bitwise OR operation.
    Pitfalls: Normalize by genome equivalence.

    #---------------------#
    `bins_x_contigs`
    7. Simply tells us which contigs got assigned to which bins.

    `bin_x_pfams`
    8. The dot product of contigs_x_bins (M7) with contigs_x_pfams (M3)

    `pfam_interesting_traits`
    9. A hand-made matrix of interesting PFAMS, description, type, pathway

    `filtered_bin_x_pfams`
    10. Take only PFAMS that have been measured and are present in the traits list.
        Add a second index so we see what is the description of the metabolism of the PFAM.
    -> Final outputted matrix for per bin metabolism.
    """

    short_name = "hit_profile"

    all_paths = """
    /norm_samples_x_pfams.tsv
    /filtered_bin_x_pfams.tsv
    """

    def __nonzero__(self): return bool(self.p.norm_samples_x_pfams)

    def __init__(self, assembly, base_dir):
        # Save parent #
        self.parent, self.assembly = assembly, assembly
        self.binner = self.assembly.results.binner
        # Auto paths #
        self.base_dir = base_dir
        self.p = AutoPaths(self.base_dir, self.all_paths)

    def run(self, verbose=False):
        # The first table #
        if verbose: print "Computing the normalized matrix samples X pfams"
        self.norm_samples_x_pfams.to_csv(self.p.norm_samples_x_pfams.path, sep='\t', float_format='%.5g')
        # The second table #
        if verbose: print "Computing the filtered matrix bins X pfams"
        self.filtered_bin_x_pfams.to_csv(self.p.filtered_bin_x_pfams.path, sep='\t', float_format='%.5g')

    @property_cached
    def samples_x_contigs(self):
        return self.binner.coverage_matrix

    @property_cached
    def contigs_x_pfams(self):
        """The results from the HMM searches for all contigs of this merged assembly."""
        # Parse #
        hmm = pandas.concat(b.pfams.hits for b in self.binner.results.bins if b)
        # Take only the top hit for each protein #
        hmm = hmm.sort_values('e_value')
        hmm = hmm.drop_duplicates('target_name')
        # Drop small e-values #
        hmm = hmm[hmm['e_value'] <= 1e-10]
        # Build new dataframe #
        result = defaultdict(lambda: defaultdict(int))
        for index, row in hmm.iterrows():
            # Remove the last number indicating proteins on the same contig #
            contig = '_'.join(row['target_name'].split('_')[:-1])
            pfam   = row['query_accession'].split('.')[0]
            result[contig][pfam] += 1
        # Fill the holes #
        result = pandas.DataFrame(result)
        result = result.fillna(0)
        result = result.astype(int)
        # Return #
        return result

    @property_cached
    def samples_x_pfams(self):
        # Drop contigs that did not have any pfams #
        df1 = self.samples_x_contigs.loc[self.contigs_x_pfams.columns]
        df2 = self.contigs_x_pfams
        # Multiply them (dot product) contigs_x_pfams and samples_x_contigs #
        assert all(df2.columns == df1.index)
        return df2.dot(df1)

    @property_cached
    def single_copy_pfams(self):
        # Load the reference file the single copy pfams #
        return pandas.io.parsers.read_csv(current_dir + 'single_copy_pfams.txt',
                                          comment   = "#",
                                          engine    = 'python',
                                          header    = None,
                                          index_col = False,
                                          squeeze   = True)
    @property_cached
    def normalization_vector(self):
        # Multiply the final matrix with this to compensate for genome equivalents across samples #
        return self.samples_x_pfams.loc[self.single_copy_pfams].median()

    @property_cached
    def norm_samples_x_pfams(self):
        # Divide by the normalization_vector #
        return self.samples_x_pfams / self.normalization_vector

    #-------------------------------------------------------------------------#
    @property_cached
    def bins_x_contigs(self):
        result = defaultdict(lambda: defaultdict(int))
        for b in self.binner.results.bins:
            for c in b.contigs:
                result[b.name][c.name] += 1
        result = pandas.DataFrame(result)
        result = result.fillna(0)
        result = result.astype(int)
        return result

    @property_cached
    def bin_x_pfams(self):
        # Drop contigs that did not have any pfams #
        df1 = self.bins_x_contigs.loc[self.contigs_x_pfams.columns]
        df2 = self.contigs_x_pfams
        # Multiply them (dot product) contigs_x_pfams and bins_x_contigs #
        assert all(df2.columns == df1.index)
        return df2.dot(df1)

    @property_cached
    def pfam_interesting_traits(self):
        # Load the file which contains special interesting pfams #
        return pandas.io.parsers.read_csv(current_dir + 'pfam_interesting_traits.csv',
                                          sep       = ",",
                                          comment   = "#",
                                          engine    = 'python')

    @property_cached
    def filtered_bin_x_pfams(self):
        # Take only the columns we like #
        pfams_with_traits = set(self.pfam_interesting_traits['PFAM'])
        pfams_we_have     = set(self.bin_x_pfams.index)
        pfams_remaining   = pfams_with_traits & pfams_we_have
        pfams_lost        = pfams_with_traits - pfams_we_have
        result            = self.bin_x_pfams.loc[pfams_remaining]
        return result

    @property_cached
    def results(self):
        results = HitProfileResults(self)
        if not results: raise Exception("You can't access results from the hit profile before running the algorithm.")
        return results

###############################################################################
class HitProfileResults(object):

    def __nonzero__(self): return bool(self.p.norm_samples_x_pfams)

    def __init__(self, hp):
        self.hp, self.parent = hp, hp
        self.p = hp.p

    @property_cached
    def norm_samples_x_pfams(self):
        return pandas.io.parsers.read_csv(self.p.norm_samples_x_pfams.path,
                                          sep='\t', index_col=0, encoding='utf-8')

    @property_cached
    def filtered_bin_x_pfams(self):
        return pandas.io.parsers.read_csv(self.p.filtered_bin_x_pfams.path,
                                          sep='\t', index_col=0, encoding='utf-8')