# Futures #
from __future__ import division

# Built-in modules #

# Internal modules #

# First party modules #
from plumbing.cache     import property_cached
from plumbing.autopaths import AutoPaths

# Third party modules #

###############################################################################
class HitProfile(object):
    """The idea is to build a profile form protein matches (such as Pfam or Tigrfam) that shows the
    different metabolic processes occurring a long a profile of samples (time series or geospatial series)

    1. Input: contigs_filtered.fasta and pfam/hits.hmmout
    2. Matrix of average coverage per contigs.   ROW: samples  - COL: Contigs - VAL: average coverage.
    3. Matrix of HMM search results.             ROW: proteins - COL: Pfams   - VAL: bit scores.
    -> From matrix 3, only take the top hit.
    4. Matrix of which contig had which protein. ROW: contigs  - COL: proteins - VAL: presence of absence.
    5: Multiply M4xM2xM3: Samples x Pfams,       ROW: samples  - COL: Pfams    - VAL: Sums of average coverages.
    -> Merge columns which have have the same Pfam by addition (Pfam with their function)
    6. Select the columns that correspond to the single copy Pfams. Calculate average by row (sample).
    -> Normalize by row using this number.

    Pitfalls: coverage variation over contigs
    Pitfalls: PFAM subfamilies (the extra number) can be thrown away and results merged by bitwise OR operation
    """

    short_name = "hit_profile"

    all_paths = """
    /lorem
    """

    def __init__(self, assembly, base_dir):
        # Save parent #
        self.parent, self.assembly = assembly, assembly
        # Auto paths #
        self.base_dir = base_dir
        self.p = AutoPaths(self.base_dir, self.all_paths)

    @property_cached
    def avgcov_x_contig(self):
        return self.assembly.results.binner.coverage_matrix

    @property_cached
    def prots_x_pfam(self):
        pass

    @property_cached
    def avgcov_x_contig(self):
        pass

    @property_cached
    def avgcov_x_contig(self):
        pass

    @property_cached
    def results(self):
        results = HitProfileResults(self)
        if not results: raise Exception("You can't access results from the hit profile before running the algorithm.")
        return results

###############################################################################
class HitProfileResults(object):
    pass