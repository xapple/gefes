# Built-in modules #

# Internal modules #
import os, inspect
from collections import defaultdict

# First party modules #
from plumbing.cache     import property_cached
from plumbing.autopaths import AutoPaths
from fasta import FASTA
from seqsearch import SeqSearch
from seqsearch.blast import BLASTdb

# Third party modules #
import pandas

# Current directory #
filename    = inspect.getframeinfo(inspect.currentframe()).filename
current_dir = os.path.dirname(os.path.abspath(filename)) + '/'

###############################################################################
class TraitAnnotations(object):
    """Let's search for a specific family of proteins in the all the FAA files
    of all the good bins."""

    short_name = "trait_annotations"
    e_value    = 1e-3
    seq_type   = "prot"

    all_paths = """
    /good_faas.fasta
    /good_faas.fasta.pin
    /stdout.txt
    /stderr.txt
    /hits.xml
    /traits_x_contigs.tsv
    /traits_x_bins.tsv
    """

    def __nonzero__(self): return bool(self.p.hits)

    def __init__(self, assembly, base_dir):
        # Save parent #
        self.parent, self.assembly = assembly, assembly
        self.binner        = self.assembly.results.binner
        self.contig_to_bin = self.binner.results.contig_id_to_bin_id
        # Auto paths #
        self.base_dir = base_dir
        self.p = AutoPaths(self.base_dir, self.all_paths)

    @property_cached
    def operons(self):
        return FASTA(current_dir + 'iron_oxidation_operons.fasta')

    @property_cached
    def good_faas(self):
        """A collection all predicted proteins in the good bins only."""
        faa = FASTA(self.p.good_fasta)
        if not faa:
            print "Collecting all proteins from good bins into a FASTA file..."
            faa.create()
            for b in self.binner.results.good_bins: faa.add(b.faa)
            faa.close()
        return faa

    @property_cached
    def good_faas_db(self):
        """A blastable database of all predicted proteins in the good bins only."""
        db = BLASTdb(self.good_faas, self.seq_type)
        if not self.p.pin:
            print "Building BLASTable database with all proteins..."
            db.makeblastdb()
        return db

    @property_cached
    def search(self):
        """The sequence similarity search to be run."""
        return SeqSearch(
              input_fasta = self.operons,
              database    = self.good_faas_db,
              seq_type    = self.seq_type,
              algorithm   = "blast",
              filtering   = {'e_value': self.e_value},
              out_path    = self.p.hits,
              _out        = self.p.stdout,
              _err        = self.p.stderr,
        )

    def run(self, verbose=True):
        # Run the search #
        self.search
        if verbose: print "Starting BLAST search."
        self.search.run()
        # Make the TSVs #
        self.traits_x_contigs.to_csv(self.p.traits_x_contigs.path, sep='\t', float_format='%.5g')
        self.traits_x_bins.to_csv(   self.p.traits_x_bins.path,    sep='\t', float_format='%.5g')
        # Make the TSV #
        return self.p.traits_x_bins

    @property_cached
    def traits_x_contigs(self):
        """A matrix with traits as columns and contigs as rows. Values represent e-values.
        For some reason BLAST adds a 'gnl|BL_ORD_ID|36683' to each title.
        Hits to the same contig from the same trait: one is picked randomly, others dropped."""
        # Build a new frame #
        result = defaultdict(lambda: defaultdict(lambda: 999999))
        for query in self.search.results:
            trait_name = query.query
            for algn, desc in zip(query.alignments, query.descriptions):
                target_name  = str(algn.hit_def.split()[0])
                contig_name  = target_name.split('_')[0]
                e_value      = desc.e
                result[trait_name][contig_name] = min(e_value, result[trait_name][contig_name])
        # Return #
        result = pandas.DataFrame(result)
        return result

    @property_cached
    def traits_x_bins(self):
        """A matrix with traits as columns and bins as rows. Values represent e-values."""
        # Build a new frame #
        result = defaultdict(lambda: defaultdict(lambda: 999999))
        for contig_name, row in self.traits_x_contigs.iterrows():
            bin_num = self.contig_to_bin[contig_name]
            for trait_name, e_value in row.iteritems():
                if pandas.isnull(e_value): continue # Get rid of NaNs
                result[trait_name][bin_num] = min(e_value, result[trait_name][bin_num])
        # Return #
        result = pandas.DataFrame(result)
        return result

    @property_cached
    def results(self):
        results = TraitAnnotationsResults(self)
        message = "You can't access results from the trait annotations before running the algorithm."
        if not results: raise Exception(message)
        return results

###############################################################################
class TraitAnnotationsResults(object):

    def __nonzero__(self): return bool(self.p.hits)

    def __init__(self, ta):
        self.ta, self.parent = ta, ta
        self.p = ta.p

    @property_cached
    def traits_x_contigs(self):
        return pandas.io.parsers.read_csv(self.p.traits_x_contigs.path,
                                          sep='\t', index_col=0, encoding='utf-8')

    @property_cached
    def traits_x_bins(self):
        return pandas.io.parsers.read_csv(self.p.traits_x_bins.path,
                                          sep='\t', index_col=0, encoding='utf-8')
