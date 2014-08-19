#!/usr/bin/env python2

"""
A script to cluster some genes.
Adapted from Moritz's undocumented "thorselia" module

The input files are on Uppmax at "/proj/b2013274/mcl".

Written by Lucas Sinclair.
Kopimi.

You can use this script from the shell like this:
$ ./gene_clustering.py
"""

# Built-in modules #
import glob
from collections import defaultdict

# Internal modules #
from gefes.common.autopaths import AutoPaths
from gefes.common.cache import property_cached
from gefes.common import natural_sort
from gefes.fasta.single import FASTA

# First party modules #
from parallelblast import BLASTdb, BLASTquery

# Third party modules #
import sh, pandas, numpy
from shell_command import shell_output

# Constants #
output_directory = "/glob/lucass/other/alex_gene_clusters"

###############################################################################
class Genes(FASTA):
    """A set of genes which are related in some way. For instance, all genes
    that clustered together when performing the MCL analysis."""

    def align(self, out_path=None):
        """We align"""
        if out_path is None: out_path = self.prefix_path + '.aln'
        sh.muscle("-in", self.path, "-out", out_path)
        #sh.gblocks(out_path, "-t=d")

    def build_tree(self):
        """We make a tree"""
        sh.raxml()

###############################################################################
class Analysis(object):
    """All sequences from all genomes are blasted against themselves.
    Then we use the MCL algorithm (Markov Cluster Algorithm) to form clusters.
    Then we count the genomes with the right number of single copy genes."""

    blast_params = {'-e': 0.1, '-W': 9, '-m': 8}

    all_paths = """
    /all_sequences.fasta
    /all_sequences.fasta.nin
    /all_sequences.blastout
    /filtered.blastout
    /filtered.abc
    /network.mci
    /dictionary.tab
    /clusters.txt
    /master_seqs.fasta
    /master_seqs.aln
    /master_seqs.tree
    """

    def __repr__(self): return '<%s object with %i genomes>' % \
        (self.__class__.__name__, len(self.genomes))

    def __init__(self, genomes, base_dir='.'):
        # Attributes #
        assert genomes
        self.genomes = genomes
        # Auto paths #
        self.base_dir = base_dir
        self.p = AutoPaths(self.base_dir, self.all_paths)

    @property
    def blast_db(self):
        """A blastable database of all genes"""
        if not self.p.all_nin.exists:
            shell_output('cat %s > %s' % (' '.join(self.genomes), self.p.all_fasta))
            BLASTdb(self.p.all_fasta).makeblastdb()
        return BLASTdb(self.p.all_fasta)

    @property
    def query(self):
        """The blast query to be executed"""
        return BLASTquery(self.blast_db, self.blast_db, self.blast_params, 'blastn', 'legacy')

    @property
    def blastout(self):
        """The blast results"""
        if not self.p.all_blastout.exists: self.query.run()
        return self.p.all_blastout

    @property
    def filtered(self):
        """We want to check the percent identify and the coverage of the hit to the query"""
        if not self.blastout.exists:
            self.blastout.copy(self.p.filtered)
        return self.p.filtered

    @property_cached
    def clusters(self):
        """A dictionary of sets. Keys are cluster names, items are sequence short ids.
        See http://bioops.info/2011/03/mcl-a-cluster-algorithm-for-graphs/"""
        if not self.p.clusters.exists:
            shell_output("cut -f 1,2,11 %s > %s" % (self.filtered, self.filtered_abc))
            sh.mcxload("-abc", self.p.filtered_abc, "--stream-mirror", "--stream-neg-log10", "-stream-tf", "ceil(200)", "-o", self.p.network, "-write-tab", self.p.dictionary)
            sh.mcl(self.p.network, "-use-tab", self.p.dictionary, "-o", self.p.clusters)
        return {"cluster_%i" % i: frozenset(line.split()) for i, line in enumerate(self.p.clusters)}

    @property_cached
    def count_table(self):
        """Genomes as columns and clusters as rows"""
        result = defaultdict(lambda: defaultdict(int))
        for genome in self.genomes:
            for c_name, c_ids in self.clusters.items():
                for gene in c_ids:
                    if gene in genome:
                        result[genome.prefix][c_name] += 1
        result = pandas.DataFrame(result)
        result = result.reindex_axis(sorted(result.index, key=natural_sort))
        result = result.fillna(0)
        return result

    @property_cached
    def single_copy_clusters(self):
        """Subset of self.clusters. Which clusters appear exactly once in each genome"""
        return {row.name: self.clusters[row.name] for i, row in self.count_table.iterrows() if numpy.all(row==1)}

    @property
    def master_seqs(self):
        """Every genome will have a 'master sequence' representing
        a concatenation of all single copy genes pertaining to it
        in a fixed order. We put this in a FASTA file."""
        if not self.p.master_fasta.exists:
            with FASTA(self.p.master_fasta) as fasta:
                for genome in self.genomes:
                    seq = ''
                    for c_name, c_ids in self.single_copy_clusters.items():
                        id_num = [i for i in c_ids if i in genome][0]
                        seq += str(genome[id_num].seq)
                    fasta.add_str(seq, name=genome.prefix)
        return Genes(self.p.master_fasta)

    @property
    def master_alignment(self):
        """The master alignment"""
        if not self.p.aln.exists: self.master_seqs.align(self.p.aln)
        return self.p.aln

    @property
    def master_tree(self):
        """The master tree"""
        if not self.p.aln.exists: self.master_seqs.align(self.p.aln)
        return self.p.aln

###############################################################################
# Real input #
files = "/proj/b2013274/mcl/*.fna"
genomes = [FASTA(path) for path in glob.glob(files)]
cluster = Analysis(genomes, output_directory)

# Test input #
test_files = output_directory + "/test/*.fna"
test_genomes = [FASTA(path) for path in glob.glob(test_files)]
test_cluster = Analysis(test_genomes, output_directory + "/test/")

# Main program #
if __name__ == '__main__':
    pass
    #cluster.make_database()
    #cluster.query.run()
    #cluster.filter_results()
    #cluster.mcl_cluster()
