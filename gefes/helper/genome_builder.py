# Built-in modules #
import os

# Internal modules #
from gefes.common.autopaths import AutoPaths
from gefes.common.slurm import nr_threads
from gefes.common import flatten
from gefes.fasta.single import FASTA
from gefes.common.cache import property_cached


# Third party modules #
import sh, numpy, matplotlib, pandas
import sklearn.mixture
from matplotlib import pyplot

###############################################################################
class GenomeBuilder(object):
    """Reassembling a genome from a bins and the reads"""

    all_paths = """
    /scaffolds.fasta
    /filtered.scaffold.fasta
    /contigs
    /graph_contigs.pdf
    /reads.fastq
    /reads.1.fastq
    /reads.2.fastq
    /bin_vs_reassembled.png
    """

    def __init__(self, parent):
        self.parent = parent
        self.bini = parent
        self.fwds=[str(p.fwd.path) for p in self.bini.parent.parent.parent]
        self.revs=[str(p.rev.path) for p in self.bini.parent.parent.parent]
        # Auto paths #
        self.base_dir = self.parent.p.rebuilt
        self.p = AutoPaths(self.base_dir, self.all_paths)

    def pull_reads(self):
        sh.bowtie2_build(self.parent.p.contigs,self.p.contigs)
        sh.bowtie2("-p", nr_threads, "-x",self.p.contigs,"-1", ",".join(self.fwds), "-2", ",".join(self.revs), "--al-conc", self.p.fastq, "-S", "/dev/null")

    def assemble_genome(self):
        spades_script = sh.Command("spades.py")
        spades_script("-o", self.base_dir,  "-1", self.p.reads_1, "-2", self.p.reads_2,  "-t", nr_threads, "--careful")

    def filter_assembly(self, prefilter = True, postfilter = False):
        original = FASTA(self.p.scaffolds)
        self.frame = pandas.DataFrame.from_dict({s.id : numpy.array(s.id.split("_"))[[3,5]] for s in original}, orient='index')
        self.frame = self.frame.astype(numpy.float)

        # filter out contigs smaller than the threshold set during the original binning (if prefilter is True)
        previous_min = numpy.array([c.length for c in self.parent.contigs]).min()
        self.classes = numpy.array([0]*len(self.frame) )
        if prefilter:
            self.classes[self.frame[0] > previous_min] = 1

        keepers = self.classes == max(set(self.classes))

        if postfilter:
            self.classes[keepers] = self.clust(numpy.log10(self.frame[keepers])) + max(set(self.classes))
            
        lens = self.frame[0].groupby(self.classes).max()
        self.frame[self.classes == lens.idxmax()]
        keepers = self.frame[self.classes == lens.idxmax()].index
        
        with FASTA(self.p.filtered) as filtered:
            filtered.add_seq([s for s in original if s.id in keepers])


    def clust(self, data):
        gmm = sklearn.mixture.GMM(n_components=2)
        gmm.fit(data[0])
        clust_x = gmm.predict(data[0])
        gmm.fit(data[1])
        clust_y = gmm.predict(data[1])
        return clust_x + 2*clust_y

    def contig_plot(self):
        fig = pyplot.figure()
        axes = fig.add_subplot(111)
        axes.set_xscale('log')
        axes.set_yscale('log')
        pyplot.scatter(self.frame[0], self.frame[1], c=self.classes, cmap=matplotlib.cm.rainbow)
        fig.savefig(self.p.graph)
        pyplot.close()

