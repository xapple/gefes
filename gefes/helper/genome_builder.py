# Built-in modules #
import os

# Internal modules #
from plumbing.autopaths import AutoPaths
from plumbing.slurm import nr_threads
from plumbing import flatten
from fasta import FASTA

# Third party modules #
import sh, numpy, matplotlib, pandas, sklearn
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

    def filter_assembly(self, cutoff=0.3):
        """Does XYZ"""
        original = FASTA(self.p.scaffold)
        data = pandas.DataFrame.from_dict({s.id : numpy.array(s.id.split("_"))[[3,5]] for s in original}, orient='index')
        data = data.astype(numpy.float)
        previous_min = numpy.array([c.length for c in self.parent.contigs]).min()
        data = data[data[0] > previous_min]
        classes = self.clust(numpy.log10(data))
        lens = data[0].groupby(classes).apply(sum)
        data[classes == lens.idxmax()]
        keepers = data[classes == lens.idxmax()].index
        with FASTA(self.p.filtered) as filtered:
            filtered.add_seq([s for s in original if s.id in keepers])

    def clust(self, data):
        gmm = sklearn.mixture.GMM(n_components=2)
        gmm.fit(data[0])
        clust_x = gmm.predict(data[0])
        gmm.fit(data[1])
        clust_y = gmm.predict(data[1])
        return clust_x + 2*clust_y

    def contig_plot(self, data, classes):
        fig = pyplot.figure()
        axes = fig.add_subplot(111)
        axes.set_xscale('log')
        axes.set_yscale('log')
        pyplot.scatter(data[0], data[1], c=classes, cmap=matplotlib.cm.rainbow)
        fig.savefig(self.p.graph)

    def graphs(self):
        return 0
