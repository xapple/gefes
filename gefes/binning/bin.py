# Futures #
from __future__ import division

# Built-in modules #

# Internal modules #
from gefes.annotation.cogs     import SingleCOGs
from gefes.annotation.checkm   import Checkm
from gefes.annotation.hmmer    import HmmQuery
from gefes.report.bin          import BinReport
from gefes.running.bin_runner  import BinRunner
from gefes.binning             import bin_graphs as graphs

# First party modules #
from plumbing.autopaths import AutoPaths
from plumbing.cache     import property_cached
from plumbing.tmpstuff  import new_temp_path
from fasta import FASTA

# Third party modules #

###############################################################################
class Bin(object):
    """A Bin is a collection of Contigs that were identified as potentially
    coming from the same population/species/strain of organisms."""

    all_paths = """
    /contigs.fasta
    /proteins.faa
    /annotation/
    /evaluation/
    /graphs/
    /report/report.pdf
    /pfam/hits.hmmout
    /tigrfam/hits.hmmout
    """

    def __str__(self):     return self.name
    def __repr__(self):    return '<%s object "%s">' % (self.__class__.__name__, self.name)
    def __len__(self):     return len(self.contig_ids)
    def __iter__(self):    return iter(self.contigs)
    def __nonzero__(self): return len(self.faa) > 0

    def __init__(self, binner, contig_ids, result_dir=None, num=None, name=None):
        """You have to specify one of either a `name` or a `num`."""
        # Save Attributes #
        self.binner     = binner
        self.contig_ids = contig_ids
        self.assembly   = self.binner.assembly
        # Num #
        if num is None: self.num = 0
        else:           self.num = num
        # Name #
        if name is None: self.name = str(self.num)
        else:            self.name = name
        # Base directory #
        if result_dir is None: self.result_dir = self.binner.p.bins_dir
        else:                  self.result_dir = result_dir
        # Auto paths #
        self.base_dir = self.result_dir + self.name + '/'
        self.p = AutoPaths(self.base_dir, self.all_paths)
        # Runner #
        self.runner = BinRunner(self)

    @property_cached
    def contigs(self):
        """A list of Contig objects."""
        return [self.assembly.results.contig_id_to_contig[c_id] for c_id in self.contig_ids]

    @property_cached
    def fasta(self):
        """A fasta file containing only the contigs pertaining to this bin in it."""
        fasta = FASTA(self.p.fasta)
        if not fasta:
            with fasta as handle:
                for contig in self.contigs:
                    handle.add_seq(contig.record)
        return fasta

    @property_cached
    def good(self):
        """Is this bin part of the <good bins> ?"""
        return self in self.binner.results.good_bins

    @property_cached
    def evaluation(self):
        """The results from evaluating the bin completeness and other metrics."""
        return Checkm(self, self.p.evaluation_dir)

    @property_cached
    def faa(self):
        """A fasta file containing only the predicted proteins
        from all contigs in this bin (through prodigal)."""
        faa = FASTA(self.p.faa)
        if not faa.exists:
            temp = FASTA(new_temp_path())
            temp.create()
            for contig in self.contigs: temp.add(contig.proteins.results.faa)
            temp.close()
            temp.move_to(faa)
        return faa

    @property_cached
    def pfams(self):
        """Use the faa file with the pfams database and hmmsearch."""
        return HmmQuery(self.faa, 'pfam', out_path=self.p.pfam_hits)

    @property_cached
    def tigrfams(self):
        """Use the faa file with the tigrfams database and hmmsearch."""
        return HmmQuery(self.faa, 'tigrfam', out_path=self.p.tigrfam_hits)

    @property_cached
    def assignment(self):
        """The imputed assignment by Phylophlan."""
        return self.binner.results.taxonomy.results.assignments.get('bin_' + self.name)

    @property_cached
    def average_coverage(self):
        """The average of coverage of this bin across all samples"""
        contig_name_to_length = lambda r: len(self.assembly.results.contig_id_to_contig[r.name])
        frame = self.binner.coverage_matrix
        frame = frame.loc[self.contig_ids]
        frame = frame.apply(lambda r: r*contig_name_to_length(r), axis=1)
        nucleotides = frame.sum().sum()
        return nucleotides / self.assembly.results.total_bp

    @property_cached
    def single_cogs(self):
        """The results from finding single copy COGs in the bin."""
        pass
        return SingleCOGs(self, self.p.annotation_dir)

    @property_cached
    def graphs(self):
        """Sorry for the black magic. The result is an object whose attributes
        are all the graphs found in bin_graphs.py initialized with this instance as
        only argument."""
        class Graphs(object): pass
        result = Graphs()
        for graph in graphs.__all__:
            cls = getattr(graphs, graph)
            setattr(result, cls.short_name, cls(self))
        return result

    @property_cached
    def report(self):
        return BinReport(self)