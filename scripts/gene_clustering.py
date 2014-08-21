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

# Futures #
from __future__ import division

# Built-in modules #
import os, glob
from collections import defaultdict

# Internal modules #
from plumbing.autopaths import AutoPaths, FilePath
from plumbing.cache import property_cached
from plumbing.common import natural_sort
from fasta import FASTA, AlignedFASTA
from parallelblast import BLASTdb, BLASTquery
from parallelblast.results import tabular_keys

# Third party modules #
import sh, pandas, numpy
from shell_command import shell_output
from Bio import Phylo
from tqdm import tqdm

# Constants #
home = os.environ['HOME'] + '/'

###############################################################################
class Cluster(object):
    """A set of genes which are related in some way. For instance, all genes
    that clustered together when performing the MCL analysis.
    (all variations of a some single copy gene class)."""

    def __repr__(self): return '<%s object number %i>' % (self.__class__.__name__, self.num)

    def __init__(self, num, ids, analysis, name=None):
        self.num = num
        self.ids = ids
        self.analysis = analysis
        self.name = "cluster_%i" % num if name is None else name
        self.path = self.analysis.p.clusters_dir + self.name + '.fasta'

    @property
    def sequences(self):
        for id_num in self.ids: yield self.analysis.blast_db[id_num]

    @property
    def fasta(self):
        """The fasta file containing the sequences of this cluster"""
        fasta = FASTA(self.path)
        if not fasta:
            with fasta as handle:
                for seq in self.sequences: handle.add_str(str(seq.seq), name=seq.id)
        return fasta

    @property
    def alignment(self):
        """The fasta file aligned with muscle"""
        alignment = AlignedFASTA(self.fasta.prefix_path + '.aln')
        if not alignment.exists:
            self.fasta.align(alignment)
            alignment.gblocks()
        return alignment

    @property
    def tree(self):
        """The tree built with raxml"""
        tree = FilePath(self.alignment.prefix_path + '.tree')
        if not tree.exists: self.alignment.build_tree(tree)
        return tree

    @property
    def phylogeny(self):
        """We can parse it with biopython"""
        return Phylo.read(self.tree.path, 'newick')

###############################################################################
class MasterCluster(Cluster):
    """Every genome will have a 'master alignment' representing
    a concatenation of all single copy genes alignments pertaining to it
    in a fixed order. We put all these alignments in the master cluster."""

    def __init__(self, analysis): self.analysis = analysis
    @property
    def sequences(self): raise NotImplementedError("MasterCluster")
    @property
    def fasta(self): raise NotImplementedError("MasterCluster")

    @property
    def alignment(self):
        alignment = AlignedFASTA(self.analysis.p.master_aln)
        if not alignment:
            with alignment as handle:
                for genome in tqdm(self.analysis.genomes):
                    seq = ''
                    for c in self.analysis.single_copy_clusters:
                        (id_num,) = (c.ids & genome.ids)
                        seq += str(c.alignment.sequences[id_num].seq)
                    handle.add_str(seq, name=genome.prefix)
        return alignment

###############################################################################
class Analysis(object):
    """All sequences from all genomes are blasted against themselves.
    Then we use the MCL algorithm (Markov Cluster Algorithm) to form clusters.
    Then we count the genomes with the right number of single copy genes."""

    blast_params = {'-e': 0.1, '-W': 9, '-m': 8}
    minimum_identity = 30.0
    mimimum_coverage = 50.0

    all_paths = """
    /all_sequences.fasta
    /all_sequences.fasta.nin
    /all_sequences.blastout
    /filtered.blastout
    /filtered.abc
    /network.mci
    /dictionary.tab
    /clusters.txt
    /master.aln
    /master.tree
    /clusters/
    """

    def __repr__(self): return '<%s object with %i genomes>' % \
        (self.__class__.__name__, len(self.genomes))

    def __init__(self, genomes, base_dir='.'):
        # Attributes #
        self.genomes = genomes
        # Auto paths #
        self.base_dir = base_dir
        self.p = AutoPaths(self.base_dir, self.all_paths)

    @property
    def blast_db(self):
        """A blastable database of all genes"""
        assert self.genomes
        if not self.p.all_nin:
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
        if not self.p.all_blastout: self.query.run()
        return self.p.all_blastout

    @property
    def filtered(self):
        """We want to check the percent identify and the coverage of the hit to the query"""
        def good_iterator(blastout):
            for line in blastout:
                info = dict(zip(tabular_keys, line.split()))
                if float(info['perc_identity']) < self.minimum_identity: continue
                query_cov = (float(info['query_end']) - float(info['query_start']))
                query_cov = 100.0 * abs(query_cov / len(self.blast_db.sql[info['query_id']]))
                subj_cov  = (float(info['subject_end']) - float(info['subject_start']))
                subj_cov  = 100.0 * abs(subj_cov  / len(self.blast_db.sql[info['subject_id']]))
                coverage = max(query_cov, subj_cov)
                if coverage < self.mimimum_coverage: continue
                yield line
        if not self.p.filtered_blastout:
            self.p.filtered_blastout.writelines(good_iterator(self.blastout))
        return self.p.filtered_blastout

    @property
    def percent_filtered(self):
        """How many hits did we filter away ?"""
        percentage = lambda x,y: (len(x)/len(y))*100 if len(y) != 0 else 0
        return "%.1f%%" % (100 - percentage(self.filtered, self.blastout))

    @property_cached
    def clusters(self):
        """A list of Clusters. See http://bioops.info/2011/03/mcl-a-cluster-algorithm-for-graphs/"""
        if not self.p.clusters.exists:
            shell_output("cut -f 1,2,11 %s > %s" % (self.filtered, self.p.filtered_abc))
            sh.mcxload("-abc", self.p.filtered_abc, "--stream-mirror", "--stream-neg-log10", "-stream-tf", "ceil(200)", "-o", self.p.network, "-write-tab", self.p.dictionary)
            sh.mcl(self.p.network, "-use-tab", self.p.dictionary, "-o", self.p.clusters)
        return [Cluster(i, frozenset(line.split()), self) for i, line in enumerate(self.p.clusters)]

    @property_cached
    def count_table(self):
        """Genomes as columns and clusters as rows"""
        result = defaultdict(lambda: defaultdict(int))
        for genome in self.genomes:
            for cluster in self.clusters:
                for id_num in cluster.ids:
                    if id_num in genome:
                        result[genome.prefix][cluster.name] += 1
        result = pandas.DataFrame(result)
        result = result.reindex_axis(sorted(result.index, key=natural_sort))
        result = result.fillna(0)
        return result

    @property_cached
    def single_copy_clusters(self):
        """Subset of self.clusters. Which clusters appear exactly once in each genome"""
        names = [row.name for i, row in self.count_table.iterrows() if numpy.all(row==1)]
        return [[c for c in self.clusters if c.name==name][0] for name in names]

    @property_cached
    def master_cluster(self): return MasterCluster(self)

###############################################################################
output_directory = home + "glob/lucass/other/alex_gene_clusters/"

# Real input #
files = "/proj/b2013274/mcl/*.fna"
genomes = [FASTA(path) for path in glob.glob(files)]
analysis = Analysis(genomes, output_directory)

# Test input #
test_files = output_directory + "/test/*.fna"
test_genomes = [FASTA(path) for path in glob.glob(test_files)]
test_analysis = Analysis(test_genomes, output_directory + "/test/")

# Main program #
#if __name__ == '__main__': Phylo.draw_ascii(analysis.master_cluster.phylogeny)