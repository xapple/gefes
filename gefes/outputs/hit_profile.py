# Futures #
from __future__ import division

# Built-in modules #
from collections import defaultdict

# Internal modules #

# First party modules #
from plumbing.cache     import property_cached
from plumbing.autopaths import AutoPaths

# Third party modules #
import pandas

###############################################################################
class HitProfile(object):
    """The idea is to build a profile form protein matches (such as Pfam or Tigrfam) that shows the
    different metabolic processes occurring a long a profile of samples (time series or geospatial series)

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

    `lorem`
    5. Select the columns that correspond to the single copy Pfams. Calculate average by row (sample).
    -> Normalize by row using this number.

    Pitfalls: Coverage variation over contigs.
    Pitfalls: PFAM subfamilies (the extra number) can be thrown away and results merged by bitwise OR operation.
    Pitfalls: Normalize by genome equivalence.
    """

    short_name = "hit_profile"

    all_paths = """
    /samples_x_pfams.tsv
    """

    def __nonzero__(self): return bool(self.p.samples_x_pfams)

    def __init__(self, assembly, base_dir):
        # Save parent #
        self.parent, self.assembly = assembly, assembly
        self.binner = self.assembly.results.binner
        # Auto paths #
        self.base_dir = base_dir
        self.p = AutoPaths(self.base_dir, self.all_paths)

    def run(self):
        content = self.samples_x_pfams.to_csv(None, sep=self.sep, float_format='%.5g')
        self.p.samples_x_pfams.writelines(content)

    @property_cached
    def samples_x_contigs(self):
        return self.binner.coverage_matrix

    @property_cached
    def contigs_x_pfams(self):
        """The results from the HMM searches for all contigs of this merged assembly"""
        # Parse #
        hmm = pandas.concat(b.pfams.hits for b in self.binner.results.bins if b)
        # Take only the top hit for each protein #
        hmm = hmm.sort_values('e_value')
        hmm = hmm.drop_duplicates('target_name')
        # Drop small e-values #
        hmm = hmm[hmm['e_value'] <= 0.001]
        # Build new dataframe #
        result = defaultdict(lambda: defaultdict(int))
        for index, row in hmm.iterrows():
            contig = row['target_name'].split('_')[0]
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
    def results(self):
        results = HitProfileResults(self)
        if not results: raise Exception("You can't access results from the hit profile before running the algorithm.")
        return results

###############################################################################
class HitProfileResults(object):

    def __nonzero__(self): return bool(self.p.samples_x_pfams)

    def __init__(self, hp):
        self.hp, self.parent = hp, hp
        self.p = hp.p

    @property_cached
    def df(self):
        return pandas.io.parsers.read_csv(self.p.samples_x_pfams.path, sep='\t', index_col=0, encoding='utf-8')