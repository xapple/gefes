#!/usr/bin/env python2

"""
A script to cluster some genes.
Adapted from Moritz's undocumented "thorsellia" module

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
from plumbing.common import natural_sort, which
from fasta import FASTA, AlignedFASTA
from parallelblast import BLASTdb, BLASTquery
from parallelblast.results import tabular_keys

# Third party modules #
import sh, pandas
from shell_command import shell_output
from Bio import Phylo
from tqdm import tqdm

# Constants #
home = os.environ['HOME'] + '/'

###############################################################################
class Genome(FASTA):
    """A FASTA file somewhere on the file system."""
    @property
    def partial(self):
        """Apparently some of them are SAGs and thus only partial."""
        return True if self.filename.startswith('2236') else False

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
    def counts(self):
        return self.analysis.count_table.loc[self.name]

    @property
    def score(self):
        """Given the genome counts, what is the single-copy likelihood score"""
        score = 0
        for name, count in self.counts.iteritems():
            partial = [g for g in self.analysis.genomes if g.prefix == name][0].partial
            if count == 0:   score +=  -5 if partial else -20
            elif count == 1: score +=  10 if partial else  10
            elif count == 2: score += -35 if partial else -30
            elif count == 3: score += -45 if partial else -40
            else: score += -100
        return score

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
        if not tree.exists:
            print "Building tree for cluster '%s'..." % self.name
            self.alignment.build_tree(tree)
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

    def __init__(self, analysis):
        self.analysis = analysis
        self.name = 'master_cluster'

    @property
    def sequences(self): raise NotImplementedError("MasterCluster")
    @property
    def fasta(self): raise NotImplementedError("MasterCluster")

    @property
    def alignment(self):
        alignment = AlignedFASTA(self.analysis.p.master_aln)
        if not alignment:
            msg = "Creating master alignment with %i genomes..."
            print msg % len(self.analysis.single_copy_clusters)
            with alignment as handle:
                for genome in self.analysis.genomes:
                    seq = ''
                    for c in tqdm(self.analysis.single_copy_clusters):
                        (id_num,) = (c.ids & genome.ids)
                        seq += str(c.alignment.sequences[id_num].seq)
                    handle.add_str(seq, name=genome.prefix)
        return alignment

###############################################################################
class Analysis(object):
    """All sequences from all genomes are blasted against themselves.
    Then we use the MCL algorithm (Markov Cluster Algorithm) to form clusters.
    Then we count the genomes with the right number of single copy genes."""

    blast_params = {'-e': 0.1, '-m': 8}
    minimum_identity = 30.0
    mimimum_coverage = 50.0
    sequence_type ='aminoacid' or 'nucleotide'

    all_paths = """
    /all_sequences.fasta
    /all_sequences.fasta.nin
    /all_sequences.fasta.pin
    /all_sequences.blastout
    /filtered.blastout
    /filtered.abc
    /network.mci
    /dictionary.tab
    /clusters.txt
    /count_table.tsv
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
        dbtype = 'nucl' if self.sequence_type == 'nucleotide' else 'prot'
        db = BLASTdb(self.p.all_fasta, dbtype=dbtype)
        if not self.p.all_nin and not self.p.all_pin:
            print "Building BLASTable database with all genes..."
            shell_output('cat %s > %s' % (' '.join(self.genomes), db))
            assert len(db.ids) == len(set(db.ids))
            db.makeblastdb()
        return db

    @property
    def query(self):
        """The blast query to be executed"""
        algorithm = 'blastn' if self.sequence_type == 'nucleotide' else 'blastp'
        return BLASTquery(self.blast_db, self.blast_db, self.blast_params, algorithm, 'legacy')

    @property
    def blastout(self):
        """The blast results"""
        if not self.p.all_blastout:
            print "Self-BLASTing database '%s'..." % self.blast_db.relative_path
            self.query.run()
        return self.p.all_blastout

    @property
    def filtered(self):
        """We want to check the percent identify and the coverage of the hit to the query"""
        def good_iterator(blastout):
            print "Filtering BLAST hits..."
            for line in tqdm(blastout, total=len(blastout)):
                info = dict(zip(tabular_keys, line.split()))
                if float(info['perc_identity']) < self.minimum_identity: continue
                query_cov = (float(info['query_end']) - float(info['query_start']))
                query_cov = 100.0 * abs(query_cov / self.blast_db.length_by_id[info['query_id']])
                subj_cov  = (float(info['subject_end']) - float(info['subject_start']))
                subj_cov  = 100.0 * abs(subj_cov  / self.blast_db.length_by_id[info['subject_id']])
                coverage = min(query_cov, subj_cov)
                if coverage < self.mimimum_coverage: continue
                yield line
        if not self.p.filtered_blastout:
            print "Making SQLite database with reads from '%s'..." % self.blastout.relative_path
            print "Result in '%s'." % self.blast_db.sql.relative_path
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
            print "Running the MCL clustering on '%s'..." % self.filtered.relative_path
            shell_output("cut -f 1,2,11 %s > %s" % (self.filtered, self.p.filtered_abc))
            sh.mcxload("-abc", self.p.filtered_abc, "--stream-mirror", "--stream-neg-log10", "-stream-tf", "ceil(200)", "-o", self.p.network, "-write-tab", self.p.dictionary)
            mcl = sh.Command(which('mcl'))
            mcl(self.p.network, "-use-tab", self.p.dictionary, "-o", self.p.clusters)
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

    def save_count_table(self):
        self.count_table = self.count_table.reindex([c.name for c in self.clusters])
        self.count_table.to_csv(str(self.p.tsv), sep='\t', encoding='utf-8')

    @property_cached
    def single_copy_clusters(self):
        """Subset of self.clusters. Which clusters appear exactly once in each genome.
        Some genomes are partial so we will be more flexible on those ones."""
        self.clusters = sorted(self.clusters, key=lambda x: x.score, reverse=True)
        return self.clusters[0:100]

    @property_cached
    def master_cluster(self): return MasterCluster(self)

###############################################################################
output_directory = home + "glob/lucass/other/alex_protein/"
extension = '.faa'

# Real input #
files = "/proj/b2013274/mcl/*" + extension
genomes = [Genome(path) for path in glob.glob(files)]
analysis = Analysis(genomes, output_directory)

# Test input #
test_files = output_directory + "/test/*" + extension
test_genomes = [Genome(path) for path in glob.glob(test_files)]
test_analysis = Analysis(test_genomes, output_directory + "/test/")

# Main program #
def run(): Phylo.draw_ascii(analysis.master_cluster.phylogeny)